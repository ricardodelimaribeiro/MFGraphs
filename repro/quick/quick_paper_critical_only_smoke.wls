SetDirectory[ParentDirectory[ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]]]];

FormatValue[x_] := Replace[x, m_Missing :> ToString[m, InputForm]];

case = "HRF Scenario 1";
d2eTiming = 0.123;
criticalTiming = 4.567;
critical = <|
  "ResultKind" -> Missing["NotAvailable"],
  "Feasibility" -> Missing["NotAvailable"],
  "Message" -> Missing["NotAvailable"],
  "KirchhoffResidual" -> Missing["NotAvailable"]
|>;

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

Export[FileNameJoin[{Directory[], "Results", "paper_critical_only_smoke.json"}], result, "RawJSON"];
Print[ExportString[result, "RawJSON", "Compact" -> False]];
