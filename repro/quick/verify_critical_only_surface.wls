repo = DirectoryName[$InputFileName, 3];
Get[FileNameJoin[{repo, "MFGraphs", "MFGraphs.wl"}]];
report = TestReport[{FileNameJoin[{repo, "MFGraphs", "Tests", "package-loading.mt"}]}];
checks = <|
    "HasCriticalCongestionSolver" -> NameQ["MFGraphs`CriticalCongestionSolver"],
    "HasSolveMFG" -> NameQ["MFGraphs`SolveMFG"],
    "HasNonLinearSolver" -> NameQ["MFGraphs`NonLinearSolver"],
    "HasMonotoneSolverFromData" -> NameQ["MFGraphs`MonotoneSolverFromData"],
    "HasTimeDependentExampleContext" -> MemberQ[$Packages, "MFGraphs`Examples`TimeDependentExamples`"],
    "HasTimeDependentSolverContext" -> MemberQ[$Packages, "MFGraphs`TimeDependentSolver`"]
|>;
failed = Length[report["FailureResults"]] + Length[report["TestsFailedWrongResultsKeys"]] + Length[report["TestsFailedWithMessagesKeys"]] + Count[
    {
        checks["HasCriticalCongestionSolver"],
        checks["HasSolveMFG"],
        !checks["HasNonLinearSolver"],
        !checks["HasMonotoneSolverFromData"],
        !checks["HasTimeDependentExampleContext"],
        !checks["HasTimeDependentSolverContext"]
    },
    False
];
Print[<|"PackageLoadingReport" -> report, "Checks" -> checks|>];
Print["CRITICAL_ONLY_SURFACE_TESTS_FAILED=", failed];
If[failed > 0, Exit[1], Exit[0]];
