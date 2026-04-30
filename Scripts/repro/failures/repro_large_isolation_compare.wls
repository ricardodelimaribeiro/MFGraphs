SetDirectory[ParentDirectory[ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]]]];
PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
Get[FileNameJoin[{Directory[], "Scripts", "BenchmarkHelpers.wls"}]];

key = "Grid0404";
data = GetExampleData[key];

Print["=== InProcess ==="];
$SafeExecuteMode = "InProcess";
inproc = SafeExecute[DataToEquations[data], 60];
Print[inproc // InputForm];

Print["=== Process ==="];
$SafeExecuteMode = "Process";
proc = SafeExecute[DataToEquations[data], 60];
Print[proc // InputForm];
