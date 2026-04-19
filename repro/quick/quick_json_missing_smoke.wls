SetDirectory[ParentDirectory[ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]]]];

FormatValue[x_] := Replace[x, m_Missing :> ToString[m, InputForm]];

result = Map[
  FormatValue,
  <|
    "Present" -> 1,
    "Absent" -> Missing["NotAvailable"],
    "Label" -> "smoke"
  |>
];

Export[FileNameJoin[{Directory[], "Results", "json_missing_smoke.json"}], result, "RawJSON"];
Print[ExportString[result, "RawJSON", "Compact" -> False]];
