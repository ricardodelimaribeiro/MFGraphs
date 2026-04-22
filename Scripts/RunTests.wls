#!/usr/bin/env wolframscript

Print["MFGraphs Test Runner"];
Print["==================="];

$RepoRoot = "/Users/ribeirrd/Documents/GitHub/MFGraphs";

PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

$TestsDir = FileNameJoin[{$RepoRoot, "MFGraphs", "Tests"}];

$TestSuites = <|
    (* Active scenario-first gates: avoid tests still tied to GetExampleData. *)
    "fast" -> {
        "graphics-public-api.mt",
        "exit-pass-through.mt",
        "exit-pass-through-nonzero.mt",
        "make-unknowns.mt"
    },
    "slow" -> {
        "run-tests-smoke.mt"
    },
    (* Legacy suites kept for manual migration tracking only. *)
    "legacy-fast" -> {
        "package-loading.mt",
        "critical-basic-cases.mt",
        "critical-numeric-backend.mt",
        "exact-mode.mt",
        "infeasible-status.mt",
        "scenario-kernel.mt",
        "solve-mfg.mt",
        "solver-contracts.mt"
    },
    "legacy-slow" -> {
        "J9F-critical_congestion.mt",
        "Jm9-critical_congestion.mt"
    }
|>;
$TestSuites["all"] = Join[$TestSuites["fast"], $TestSuites["slow"]];
$TestSuites["legacy"] = Join[$TestSuites["legacy-fast"], $TestSuites["legacy-slow"]];
$TestSuites["full"] = Join[$TestSuites["all"], $TestSuites["legacy"]];

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
Print["Files: ", Length[suiteFiles]];
Print[];

$totalPassed = 0;
$totalFailed = 0;
$allFailed = {};

Do[
    Print["--- ", file, " ---"];
    $t0 = AbsoluteTime[];
    report = TestReport[FileNameJoin[{$TestsDir, file}]];
    $elapsed = Round[AbsoluteTime[] - $t0];
    passed = report["TestsSucceededCount"];
    failed = report["TestsFailedCount"];
    $totalPassed += passed;
    $totalFailed += failed;
    Print["  Passed: ", passed, "  Failed: ", failed, "  (", $elapsed, "s)"];
    If[failed > 0,
        failedHere = Select[report["TestResults"], !MatchQ[#["Outcome"], "Success"] &];
        KeyValueMap[
            (Print["  [FAIL] ", #2["TestID"]]; AppendTo[$allFailed, #2["TestID"]]) &,
            failedHere
        ]
    ];
    Print[];
    ,
    {file, suiteFiles}
];

Print["==================="];
Print["Total passed: ", $totalPassed];
Print["Total failed: ", $totalFailed];

If[$totalFailed > 0,
    Print["\nFailed tests:"];
    Scan[Print["  ", #] &, $allFailed];
    Exit[1]
,
    Exit[0]
];
