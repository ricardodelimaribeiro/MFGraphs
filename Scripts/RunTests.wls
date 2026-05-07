#!/usr/bin/env wolframscript

Print["MFGraphs Test Runner"];
Print["==================="];

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];

PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

$TestsDir = FileNameJoin[{$RepoRoot, "MFGraphs", "Tests"}];

$TestSuites = <|
    (* Active suite for the scenario-kernel phase. *)
    "fast" -> {
        "scenario-kernel.mt",
        "symbolic-unknowns.mt",
        "reduce-system.mt",
        "scenario-consistency.mt",
        "graphicsTools.mt",
        "orchestration.mt",
        "dnf-reducer.mt",
        "tawaf.mt"
    },
    (* Archived tests are intentionally excluded from active CI.
       Use this suite only for explicit compatibility checks. *)
    "archive" -> {
        "archive/package-loading.mt",
        "archive/graphics-public-api.mt",
        "archive/critical-basic-cases.mt",
        "archive/critical-numeric-backend.mt",
        "archive/exit-pass-through.mt",
        "archive/exit-pass-through-nonzero.mt",
        "archive/exact-mode.mt",
        "archive/infeasible-status.mt",
        "archive/solve-mfg.mt",
        "archive/solver-contracts.mt",
        "archive/J9F-critical_congestion.mt",
        "archive/Jm9-critical_congestion.mt",
        "archive/run-tests-smoke.mt",
        "archive/inconsistent-switching-critical.mt",
        "archive/symbolic-underdetermined.mt"
    }
|>;
$TestSuites["all"] = $TestSuites["fast"];
$TestSuites["full"] = Join[$TestSuites["fast"], $TestSuites["archive"]];

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
