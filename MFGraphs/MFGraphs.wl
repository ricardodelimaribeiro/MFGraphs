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

(* Load subpackages before BeginPackage so each one establishes MFGraphs` independently. *)
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

MFGPrint::usage =
"MFGPrint[args___] prints args only when $MFGraphsVerbose is True.";

MFGPrintTemporary::usage =
"MFGPrintTemporary[args___] prints a temporary message only when $MFGraphsVerbose is True.";

GridScenario::usage =
"GridScenario[dims, entries, exits] creates a scenario on a directed GridGraph[dims]. \
{n} gives a chain with vertices 1..n; {r,c} gives an r\[Times]c grid with vertices 1..r*c (row-major). \
Optional: sc (switching costs, default {}), alpha, V, g (Hamiltonian defaults from $DefaultHamiltonian).";

CycleScenario::usage =
"CycleScenario[n, entries, exits] creates a scenario on a directed n-cycle (1->2->...->n->1), \
vertices 1..n. Optional: sc, alpha, V, g.";

GraphScenario::usage =
"GraphScenario[graph, entries, exits] creates a scenario from any WL directed Graph object. \
Optional: sc, alpha, V, g.";

AMScenario::usage =
"AMScenario[vl, am, entries, exits] creates a scenario from an explicit vertices list vl \
and adjacency matrix am. Optional: sc, alpha, V, g.";

GetExampleScenario::usage =
"GetExampleScenario[n] returns a 6-arg factory Function[{entries,exits,sc,alpha,V,g}, scenario[...]] \
for built-in example n. Topology is baked in; all parameters are caller-supplied. \
GetExampleScenario[n, entries, exits] calls the factory using the canonical switching costs \
for that case (sc=Automatic resolves via $CaseDefaultSC, defaulting to {} if none defined) \
and standard Hamiltonian defaults (alpha=1, V=0, g=Function[z,-1/z]). \
Additional optional arguments override each default in order: sc, alpha, V, g. \
Pass sc={} explicitly to force no switching costs. \
entries={{vertex,flow},...}, exits={{vertex,cost},...}, sc={{i,k,j,cost},...}. \
Returns $Failed for unknown keys.";

(* --- Private Section --- *)
Begin["`Private`"];

(* Defaults *)
$MFGraphsVerbose = False;

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
MFGPrint[args___] := If[$MFGraphsVerbose, Print[args]];
MFGPrintTemporary[args___] := If[$MFGraphsVerbose, PrintTemporary[args], Null];
(* :!CodeAnalysis::EndBlock:: *)

NumberVectorQ[j_] := VectorQ[j, NumericQ];

IsFeasible[assoc_Association] :=
    Lookup[assoc, "Feasibility", "Infeasible"] === "Feasible";

alpha[_] := 1;
Cost[m_, edge_] := m^alpha[edge];

CheckFlowFeasibility[Null] := "Infeasible";
CheckFlowFeasibility[assoc_Association] :=
    Module[{flowKeys, flowVals},
        flowKeys = Select[Keys[assoc], MatchQ[#, _j] &];
        flowVals = Lookup[assoc, flowKeys];
        If[flowVals === {} || Min[Select[flowVals, NumericQ]] < 0,
            "Infeasible", "Feasible"]
    ];

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
