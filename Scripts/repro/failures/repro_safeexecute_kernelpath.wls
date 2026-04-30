SetDirectory[ParentDirectory[ParentDirectory[DirectoryName[ExpandFileName[$InputFileName]]]]];
PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
Get[FileNameJoin[{Directory[], "Scripts", "BenchmarkHelpers.wls"}]];

Print["=== InProcess 1+1 ==="];
$SafeExecuteMode = "InProcess";
Print[SafeExecute[1+1, 10] // InputForm];

Print["=== Process 1+1 default kernelPath ==="];
$SafeExecuteMode = "Process";
Print[SafeExecute[1+1, 10] // InputForm];

Print["=== Direct isolated 1+1 explicit kernelPath ==="];
(* Path is hardcoded intentionally: this block tests the explicit override mechanism,
   not the auto-resolution path. Adjust to your local wolframscript location if needed. *)
Print[SafeExecuteIsolated[1+1, 10, "KernelPath" -> "/usr/local/bin/wolframscript"] // InputForm];
