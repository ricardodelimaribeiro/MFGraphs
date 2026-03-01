(* Wolfram Language Package *)

(* Created by the Wolfram Workbench May 5, 2020 *)
(* To distribute, use ideas from https://community.wolfram.com/groups/-/m/t/214901 *)

BeginPackage["MFGraphs`"];

(* Verbose flag: set to False to suppress progress messages *)
$MFGraphsVerbose::usage =
"$MFGraphsVerbose controls whether progress and timing messages are printed.
Set $MFGraphsVerbose = False to suppress output. Default is True.";
$MFGraphsVerbose = True;

MFGPrint::usage =
"MFGPrint[args___] prints args only when $MFGraphsVerbose is True.";
MFGPrint[args___] := If[$MFGraphsVerbose, Print[args]];

MFGPrintTemporary::usage =
"MFGPrintTemporary[args___] prints a temporary message only when $MFGraphsVerbose is True.";
MFGPrintTemporary[args___] := If[$MFGraphsVerbose, PrintTemporary[args], Null];

(* Load submodules in dependency order *)
Get["MFGraphs`Examples`ExamplesData`"];
Get["MFGraphs`ZAnd`"];
Get["MFGraphs`D2E2`"];
Get["MFGraphs`NonLinearSolver`"];
Get["MFGraphs`Monotone`"];

EndPackage[];
