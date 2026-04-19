SetDirectory[ParentDirectory[ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]]]];
PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];

FormatValue[x_] := Replace[x, m_Missing :> ToString[m, InputForm]];

$MFGraphsVerbose = False;
$HistoryLength = 0;
Quiet @ Check[MFGraphs`Private`ClearSolveCache[], Null];
Quiet @ Check[ClearSystemCache[], Null];

case = "HRF Scenario 1";
data = GetExampleData[case];

d2eTiming = AbsoluteTiming[d2e = DataToEquations[data];][[1]];
criticalTiming = AbsoluteTiming[critical = CriticalCongestionSolver[d2e];][[1]];

result = Map[
  FormatValue,
  <|
    "Case" -> case,
    "Tier" -> "paper",
    "IncludedSolvers" -> {"DataToEquations", "CriticalCongestion"},
    "D2ETimeSeconds" -> d2eTiming,
    "CriticalTimeSeconds" -> criticalTiming,
    "CriticalResultKind" -> Lookup[critical, "ResultKind", Missing["NotAvailable"]],
    "CriticalFeasibility" -> Lookup[critical, "Feasibility", Missing["NotAvailable"]],
    "CriticalMessage" -> Lookup[critical, "Message", Missing["NotAvailable"]],
    "KirchhoffResidual" -> Lookup[critical, "KirchhoffResidual", Missing["NotAvailable"]],
    "Timestamp" -> DateString[]
  |>
];

Export[FileNameJoin[{Directory[], "Results", "paper_critical_only.json"}], result, "RawJSON"];
Print[ExportString[result, "RawJSON", "Compact" -> False]];
