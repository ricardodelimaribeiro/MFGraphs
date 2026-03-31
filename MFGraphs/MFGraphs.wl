(* Wolfram Language Package *)
(*
   MFGraphs: A Wolfram Language package for Mean Field Games on Networks.
   
   Main Components:
   - DataToEquations: Symbolic graph-to-equation converter.
   - NonLinearSolver: Iterative solver for general congestion.
   - Monotone: Hessian Riemannian Flow numerical solver.
   - DNFReduce: Symbolic logical reduction engine.
*)

(* Created by the Wolfram Workbench May 5, 2020 *)
(* To distribute, use ideas from https://community.wolfram.com/groups/-/m/t/214901 *)

BeginPackage["MFGraphs`"];

(* Verbose flag: set to False to suppress progress messages *)
$MFGraphsVerbose::usage =
"$MFGraphsVerbose controls whether progress and timing messages are printed.
Set $MFGraphsVerbose = False to suppress output. Default is True.";
$MFGraphsVerbose = False;

$MFGraphsParallelThreshold::usage =
"$MFGraphsParallelThreshold controls the minimum list length before ParallelMap or
ParallelTable is used instead of Map or Table. Default is 6. Set to Infinity to
disable all parallelism.";
$MFGraphsParallelThreshold = 6;

$MFGraphsParallelReady::usage =
"$MFGraphsParallelReady is True after ParallelNeeds[\"MFGraphs`\"] has been called on
all subkernels. Reset to False (e.g. after LaunchKernels[]) to force re-initialization.";
$MFGraphsParallelReady = False;

MFGPrint::usage =
"MFGPrint[args___] prints args only when $MFGraphsVerbose is True.";
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
MFGPrint[args___] := If[$MFGraphsVerbose, Print[args]];

MFGPrintTemporary::usage =
"MFGPrintTemporary[args___] prints a temporary message only when $MFGraphsVerbose is True.";
MFGPrintTemporary[args___] := If[$MFGraphsVerbose, PrintTemporary[args], Null];
(* :!CodeAnalysis::EndBlock:: *)

EnsureParallelKernels::usage =
"EnsureParallelKernels[] launches parallel subkernels if none are running.";
EnsureParallelKernels[] := If[$KernelCount === 0, LaunchKernels[]];

MFGParallelMap::usage =
"MFGParallelMap[f, list] applies f to each element of list, using ParallelMap
when Length[list] >= $MFGraphsParallelThreshold, otherwise Map.";
MFGParallelMap[f_, list_List] :=
    If[Length[list] >= $MFGraphsParallelThreshold,
        (EnsureParallelKernels[]; ParallelMap[f, list]),
        f /@ list
    ];

Begin["`Private`"];

CheckFlowFeasibility[Null] := "Infeasible";
CheckFlowFeasibility[assoc_Association] :=
    Module[{flowKeys, flowVals},
        flowKeys = Select[Keys[assoc], MatchQ[#, _j] &];
        flowVals = Lookup[assoc, flowKeys];
        If[flowVals === {} || Min[Select[flowVals, NumericQ]] < 0,
            "Infeasible", "Feasible"]
    ];

SelectFlowAssociation[assoc_Association] :=
    Association @ KeySelect[assoc, MatchQ[#, _j] &];

BuildSolverComparisonData[Eqs_Association, solution_] :=
    Module[{missing, flowAssoc, edgeList, edgePairs, signedFlowRules, signedExprs,
         signedVals, B, KM, jj, jjVals, residual},
        missing = Missing["NotAvailable"];
        edgeList = Lookup[Eqs, "edgeList", {}];
        If[!AssociationQ[solution],
            Return[
                <|
                    "ComparableEdges" -> edgeList,
                    "FlowAssociation" -> missing,
                    "SignedEdgeFlows" -> missing,
                    "ComparableFlowVector" -> missing,
                    "KirchhoffResidual" -> missing
                |>,
                Module
            ]
        ];
        flowAssoc = SelectFlowAssociation[solution];
        edgePairs = List @@@ edgeList;
        signedFlowRules = Lookup[Eqs, "SignedFlows", <||>];
        signedExprs =
            (If[KeyExistsQ[signedFlowRules, #], signedFlowRules[#], missing] & /@ edgePairs) /.
                Join[
                    Lookup[Eqs, "RuleBalanceGatheringFlows", <||>],
                    Lookup[Eqs, "RuleExitFlowsIn", <||>],
                    Lookup[Eqs, "RuleEntryOut", <||>]
                ];
        signedVals = Quiet @ Check[signedExprs /. flowAssoc, missing];
        If[!ListQ[signedVals] || !VectorQ[signedVals, NumericQ],
            signedVals = missing
        ];
        {B, KM, jj} = MFGraphs`GetKirchhoffLinearSystem[Eqs];
        jjVals = Lookup[flowAssoc, jj, missing];
        residual =
            If[ListQ[jjVals] && VectorQ[jjVals, NumericQ],
                N @ Norm[KM . jjVals - B, Infinity],
                missing
            ];
        <|
            "ComparableEdges" -> edgeList,
            "FlowAssociation" -> flowAssoc,
            "SignedEdgeFlows" ->
                If[signedVals === missing, missing, AssociationThread[edgeList, signedVals]],
            "ComparableFlowVector" -> signedVals,
            "KirchhoffResidual" -> residual
        |>
    ];

MakeSolverResult[
    solver_String,
    resultKind_String,
    feasibility_,
    message_,
    solution_,
    extraAssoc_: <||>
] :=
    Join[
        <|
            "Solver" -> solver,
            "ResultKind" -> resultKind,
            "Feasibility" -> feasibility,
            "Message" -> message,
            "Solution" -> solution
        |>,
        extraAssoc
    ];

End[];

(* Load submodules in dependency order *)
Get["MFGraphs`Examples`ExamplesData`"];
Get["MFGraphs`Examples`TimeDependentExamples`"];
Get["MFGraphs`DNFReduce`"];
Get["MFGraphs`DataToEquations`"];
Get["MFGraphs`NonLinearSolver`"];
Get["MFGraphs`Monotone`"];
Get["MFGraphs`TimeDependentSolver`"];

EndPackage[];
