repo = DirectoryName[$InputFileName, 3];
Get[FileNameJoin[{repo, "MFGraphs", "MFGraphs.wl"}]];
report = TestReport[{FileNameJoin[{repo, "MFGraphs", "Tests", "inconsistent-switching-critical.mt"}]}];
failed = Length[report["FailureResults"]] + Length[report["TestsFailedWrongResultsKeys"]] + Length[report["TestsFailedWithMessagesKeys"]];
Print[report];
Print["INCONSISTENT_SWITCHING_CRITICAL_FAILED=", failed];
If[failed > 0, Exit[1], Exit[0]];
