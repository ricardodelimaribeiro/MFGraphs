(* Wolfram Language Package *)
(*
   MFGraphs: A Wolfram Language package for Mean Field Games on Networks.

   Current phase: core scenario kernels.

   This top-level file is a thin bootstrap. It owns the BeginPackage/
   EndPackage block, shared verbosity/parallel helpers, and loads the
   active typed-kernel submodules listed below.

   Main Components:
   - Scenario              - typed scenario kernel.
   - Examples/ExampleScenarios
                           - named scenario constructors/factories.
   - Unknowns              - symbolic unknown-family construction.
   - System                - structural equation-system kernel.
*)

(* Created by the Wolfram Workbench May 5, 2020 *)
(* To distribute, use ideas from https://community.wolfram.com/groups/-/m/t/214901 *)

(* Scenario.wl and ExampleScenarios.wl are loaded via explicit paths before BeginPackage
   so each can establish MFGraphs` independently. Unknowns and System are loaded via
   context names inside Private after BeginPackage, which is the standard submodule pattern. *)
With[{$mfgDir = DirectoryName[$InputFileName]},
    Get[FileNameJoin[{$mfgDir, "Scenario.wl"}]];
    Get[FileNameJoin[{$mfgDir, "Examples", "ExampleScenarios.wl"}]];
    Null
];

BeginPackage["MFGraphs`"];

(* --- Public API declarations (Usage strings) --- *)

alpha::usage = "alpha[edge] is the congestion exponent for an edge. Default is 1.";
Cost::usage = "Cost[m, edge] is the congestion cost function.";

j::usage = "j[v, e] or j[v, e1, e2] represents a flow variable.";
u::usage = "u[v, e] represents a utility/potential variable.";
z::usage = "z[v] represents a vertex potential variable.";

(* Verbose flag: set to False to suppress progress messages *)
$MFGraphsVerbose::usage =
"$MFGraphsVerbose controls whether progress and timing messages are printed.
Set $MFGraphsVerbose = False to suppress output. Default is True.";

$MFGraphsParallelThreshold::usage =
"$MFGraphsParallelThreshold controls the minimum list length before ParallelMap or
ParallelTable is used instead of Map or Table. Default is 6. Set to Infinity to
disable all parallelism.";

$MFGraphsParallelReady::usage =
"$MFGraphsParallelReady is True after ParallelNeeds[\"MFGraphs`\"] has been called on
all subkernels. Reset to False (e.g. after LaunchKernels[]) to force re-initialization.";

MFGPrint::usage =
"MFGPrint[args___] prints args only when $MFGraphsVerbose is True.";

MFGPrintTemporary::usage =
"MFGPrintTemporary[args___] prints a temporary message only when $MFGraphsVerbose is True.";

EnsureParallelKernels::usage =
"EnsureParallelKernels[] launches parallel subkernels if none are running.";


(* --- Private Section --- *)
Begin["`Private`"];

(* Defaults *)
$MFGraphsVerbose = False;
$MFGraphsParallelThreshold = 6;
$MFGraphsParallelReady = False;

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
MFGPrint[args___] := If[$MFGraphsVerbose, Print[args]];
MFGPrintTemporary[args___] := If[$MFGraphsVerbose, PrintTemporary[args], Null];
(* :!CodeAnalysis::EndBlock:: *)

EnsureParallelKernels[] := If[$KernelCount === 0, LaunchKernels[]];

alpha[_] := 1;
Cost[m_, edge_] := m^alpha[edge];

(* Load submodules in dependency order *)
Scan[
    Get,
    {
        "MFGraphs`Unknowns`",
        "MFGraphs`System`"
    }
];

End[];

EndPackage[];
