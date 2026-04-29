(* Wolfram Language package *)
(* primitives.wl — shared symbolic primitives and runtime flags for MFGraphs.

   This is the base layer. All other sub-packages depend on this one.
   Mirrors the role of utilities.m in the numerics package.

   Provides:
     - Symbolic heads: j, u, z, alpha, Cost
     - Runtime flags: $MFGraphsVerbose, $MFGraphsParallelThreshold, $MFGraphsParallelReady
     - Print utilities: mfgPrint, mfgPrintTemporary, ensureParallelKernels
*)

BeginPackage["primitives`"]

j::usage = "j[a, b] is the flow on edge {a,b} from a to b. j[r, i, w] is the fraction of the flow j[r,i] that transitions to edge e_{i,w} at junction i.";
u::usage = "u[a, b] is the value of the value function at vertex b of edge {a,b}.";
z::usage = "z[v] represents a vertex potential variable.";
alpha::usage = "alpha[edge] is the congestion exponent for an edge. Default is 1.";
Cost::usage = "Cost[m, edge] is the congestion cost function.";

$MFGraphsVerbose::usage =
"$MFGraphsVerbose controls whether progress and timing messages are printed. \
Set $MFGraphsVerbose = False to suppress output. Default is False.";

$MFGraphsParallelThreshold::usage =
"$MFGraphsParallelThreshold controls the minimum list length before ParallelMap or \
ParallelTable is used instead of Map or Table. Default is 6. Set to Infinity to \
disable all parallelism.";

$MFGraphsParallelReady::usage =
"$MFGraphsParallelReady is True after ParallelNeeds[\"MFGraphs`\"] has been called on \
all subkernels. Reset to False (e.g. after LaunchKernels[]) to force re-initialization.";

mfgPrint::usage = "mfgPrint[args___] prints args only when $MFGraphsVerbose is True.";

mfgPrintTemporary::usage = "mfgPrintTemporary[args___] prints a temporary message only when $MFGraphsVerbose is True.";

ensureParallelKernels::usage = "ensureParallelKernels[] launches parallel subkernels if none are running.";

Begin["`Private`"]

$MFGraphsVerbose          = False;
$MFGraphsParallelThreshold = 6;
$MFGraphsParallelReady    = False;

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
mfgPrint[args___] := If[$MFGraphsVerbose, Print[args]];
mfgPrintTemporary[args___] := If[$MFGraphsVerbose, PrintTemporary[args], Null];
(* :!CodeAnalysis::EndBlock:: *)

ensureParallelKernels[] := If[$KernelCount === 0, LaunchKernels[]];

alpha[_] := 1;
Cost[m_, edge_] := m^alpha[edge];

End[]

EndPackage[]
