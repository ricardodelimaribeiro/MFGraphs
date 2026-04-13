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
when Length[list] >= $MFGraphsParallelThreshold (and kernels are already running)
or Length[list] >= $MFGraphsParallelLaunchThreshold (when kernels must be launched).
This avoids paying ~3s kernel launch overhead for small workloads.";

SolveMFG::usage =
"SolveMFG[input, Method -> m, opts] dispatches to the selected stationary solver and \
returns its standardized result association.

Supported methods:
  \"Automatic\" — cascades CriticalCongestion -> Monotone -> NonLinear, falling back on \
infeasibility or failure. When \"CongestionExponentFunction\" specifies alpha != 1, \
CriticalCongestion is skipped (it only handles alpha=1).
  \"CriticalCongestion\"
  \"NonLinear\"
  \"Monotone\"

Hamiltonian options (forwarded to NonLinear and Monotone solvers):
  \"PotentialFunction\" (default Automatic)
  \"CongestionExponentFunction\" (default Automatic, meaning alpha=1). Accepts a scalar, \
a Function[edge, ...], or an Association mapping edges to exponents (unlisted edges \
default to 1).
  \"InteractionFunction\" (default Automatic)

input can be either raw data (association accepted by DataToEquations) or an already \
compiled equation association.";

$MFGraphsParallelLaunchThreshold::usage =
"$MFGraphsParallelLaunchThreshold is the minimum list length required to justify
launching parallel subkernels. Only applies when $KernelCount === 0. Default is 50.";
$MFGraphsParallelLaunchThreshold = 50;

MFGParallelMap[f_, list_List] :=
    Module[{threshold},
        threshold = If[$KernelCount === 0,
            $MFGraphsParallelLaunchThreshold,
            $MFGraphsParallelThreshold
        ];
        If[Length[list] >= threshold,
            (EnsureParallelKernels[]; ParallelMap[f, list]),
            f /@ list
        ]
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

BuildBoundaryMassData[Eqs_Association, flowAssoc_Association] :=
    Module[{missing, entryPairs, exitPairs, entryFlows, exitFlows, inTotal, outTotal, residual},
        missing = Missing["NotAvailable"];
        entryPairs = List @@@ Lookup[Eqs, "entryEdges", {}];
        exitPairs = List @@@ Lookup[Eqs, "exitEdges", {}];
        entryFlows =
            If[entryPairs === {},
                {},
                Lookup[flowAssoc, j @@@ entryPairs, missing]
            ];
        exitFlows =
            If[exitPairs === {},
                {},
                Lookup[flowAssoc, j @@@ exitPairs, missing]
            ];
        inTotal =
            If[
                ListQ[entryFlows] && VectorQ[entryFlows, NumericQ],
                N @ Total[entryFlows],
                missing
            ];
        outTotal =
            If[
                ListQ[exitFlows] && VectorQ[exitFlows, NumericQ],
                N @ Total[exitFlows],
                missing
            ];
        residual =
            If[
                NumericQ[inTotal] && NumericQ[outTotal],
                N @ Abs[inTotal - outTotal],
                missing
            ];
        <|
            "BoundaryMassIn" -> inTotal,
            "BoundaryMassOut" -> outTotal,
            "BoundaryMassResidual" -> residual
        |>
    ];

BuildBoundaryMassData[_, _] :=
    <|
        "BoundaryMassIn" -> Missing["NotAvailable"],
        "BoundaryMassOut" -> Missing["NotAvailable"],
        "BoundaryMassResidual" -> Missing["NotAvailable"]
    |>;

BuildUtilityReductionResidualData[Eqs_Association, solution_Association] :=
    Module[{missing, utilityReduction, classes, classResiduals, overallResidual,
      classCount, reducedVarCount},
        missing = Missing["NotAvailable"];
        utilityReduction = Lookup[Eqs, "UtilityReduction", Missing["NotAvailable"]];
        If[!AssociationQ[utilityReduction],
            Return[
                <|
                    "UtilityClassResidual" -> missing,
                    "UtilityClassCount" -> missing,
                    "UtilityReducedVariableCount" -> missing
                |>,
                Module
            ]
        ];
        classes = Lookup[utilityReduction, "UtilityClasses", {}];
        classCount =
            If[ListQ[classes], Length[classes], missing];
        reducedVarCount = Lookup[utilityReduction, "ReducedVariableCount", missing];
        classResiduals = Cases[
            Map[
                Function[class,
                    Module[{vals},
                        vals = Quiet @ Check[N @ (class /. solution), $Failed];
                        If[ListQ[vals] && VectorQ[vals, NumericQ] && Length[vals] > 1,
                            N @ Max @ Abs[vals - First[vals]],
                            If[ListQ[vals] && VectorQ[vals, NumericQ], 0., Missing["NotAvailable"]]
                        ]
                    ]
                ],
                classes
            ],
            _?NumericQ
        ];
        overallResidual =
            If[classResiduals === {}, missing, N @ Max[classResiduals]];
        <|
            "UtilityClassResidual" -> overallResidual,
            "UtilityClassCount" -> classCount,
            "UtilityReducedVariableCount" -> reducedVarCount
        |>
    ];

BuildUtilityReductionResidualData[_, _] :=
    <|
        "UtilityClassResidual" -> Missing["NotAvailable"],
        "UtilityClassCount" -> Missing["NotAvailable"],
        "UtilityReducedVariableCount" -> Missing["NotAvailable"]
    |>;

BuildSolverComparisonData[Eqs_Association, solution_] :=
    Module[{missing, flowAssoc, edgeList, edgePairs, signedFlowRules, signedExprs,
         signedVals, B, KM, jj, jjVals, residual, numericState, encodedFlowVec,
         fastSignedVals, fastResidual, usedFastSignedQ, usedFastResidualQ,
         boundaryMassData, utilityReductionResidualData},
        missing = Missing["NotAvailable"];
        edgeList = Lookup[Eqs, "edgeList", {}];
        If[!AssociationQ[solution],
            Return[
                <|
                    "ComparableEdges" -> edgeList,
                    "FlowAssociation" -> missing,
                    "SignedEdgeFlows" -> missing,
                    "ComparableFlowVector" -> missing,
                    "KirchhoffResidual" -> missing,
                    "BoundaryMassIn" -> missing,
                    "BoundaryMassOut" -> missing,
                    "BoundaryMassResidual" -> missing,
                    "UtilityClassResidual" -> missing,
                    "UtilityClassCount" -> missing,
                    "UtilityReducedVariableCount" -> missing
                |>,
                Module
            ]
        ];
        flowAssoc = SelectFlowAssociation[solution];
        boundaryMassData = BuildBoundaryMassData[Eqs, flowAssoc];
        utilityReductionResidualData = BuildUtilityReductionResidualData[Eqs, solution];
        numericState = Lookup[Eqs, "NumericState", Missing["NotAvailable"]];
        usedFastSignedQ = False;
        usedFastResidualQ = False;
        signedVals = missing;
        residual = missing;

        If[AssociationQ[numericState],
            encodedFlowVec = Quiet @ Check[
                EncodeFlowAssociation[numericState, flowAssoc],
                missing
            ];
            If[ListQ[encodedFlowVec] && VectorQ[encodedFlowVec, NumericQ],
                fastSignedVals = Quiet @ Check[
                    ComputeSignedEdgeFlowsFast[numericState, encodedFlowVec],
                    missing
                ];
                If[
                    ListQ[fastSignedVals] &&
                    VectorQ[fastSignedVals, NumericQ] &&
                    Length[fastSignedVals] === Length[edgeList],
                    signedVals = fastSignedVals;
                    usedFastSignedQ = True
                ];
                fastResidual = Quiet @ Check[
                    ComputeKirchhoffResidualFast[numericState, encodedFlowVec],
                    missing
                ];
                If[NumericQ[fastResidual],
                    residual = fastResidual;
                    usedFastResidualQ = True
                ]
            ]
        ];

        If[!usedFastSignedQ,
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
            ]
        ];

        If[!usedFastResidualQ,
            {B, KM, jj} = MFGraphs`GetKirchhoffLinearSystem[Eqs];
            jjVals = Lookup[flowAssoc, jj, missing];
            residual =
                If[ListQ[jjVals] && VectorQ[jjVals, NumericQ],
                    N @ Norm[KM . jjVals - B, Infinity],
                    missing
                ]
        ];
        Join[
            <|
                "ComparableEdges" -> edgeList,
                "FlowAssociation" -> flowAssoc,
                "SignedEdgeFlows" ->
                    If[signedVals === missing, missing, AssociationThread[edgeList, signedVals]],
                "ComparableFlowVector" -> signedVals,
                "KirchhoffResidual" -> residual
            |>,
            boundaryMassData,
            utilityReductionResidualData
        ]
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
        If[
            AssociationQ[extraAssoc] && !KeyExistsQ[extraAssoc, "Convergence"],
            <|"Convergence" -> Missing["NotApplicable"]|>,
            <||>
        ],
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

Options[SolveMFG] = DeleteDuplicatesBy[
    Join[
        {Method -> "Automatic"},
        Options[NonLinearSolver],
        Options[MonotoneSolverFromData],
        Options[CriticalCongestionSolver]
    ],
    First
];

SolveMFGUnknownMethodResult[method_] :=
    <|
        "Solver" -> "SolveMFG",
        "ResultKind" -> "Failure",
        "Feasibility" -> Missing["NotApplicable"],
        "Message" -> "UnknownMethod",
        "Method" -> method,
        "Solution" -> Missing["NotAvailable"],
        "Convergence" -> Missing["NotApplicable"]
    |>;

SolveMFGCompiledInputQ[input_Association] :=
    KeyExistsQ[input, "EqGeneral"] &&
    KeyExistsQ[input, "js"] &&
    KeyExistsQ[input, "jts"];

$SolveMFGRequiredEnvelopeKeys = {"Solver", "ResultKind", "Feasibility", "Message"};

SolveMFGClassifyOutcome[result_] :=
    Module[{resultKind, feasibility, message},
        If[!AssociationQ[result], Return["Abort"]];
        If[!And @@ (KeyExistsQ[result, #] & /@ $SolveMFGRequiredEnvelopeKeys), Return["Abort"]];
        resultKind = Lookup[result, "ResultKind", Missing["NotAvailable"]];
        feasibility = Lookup[result, "Feasibility", Missing["NotAvailable"]];
        message = Lookup[result, "Message", Missing["NotAvailable"]];
        Which[
            resultKind === "Success" && feasibility === "Feasible", "Done",
            resultKind === "Degenerate" && message === "DegenerateCase", "Done",
            feasibility === "Infeasible" || MemberQ[{"Failure", "NonConverged"}, resultKind], "Fallback",
            True, "Abort"
        ]
    ];

SolveMFGBuildTraceEntry[method_, result_, decision_] :=
    Module[{solver, resultKind, feasibility, message},
        solver = If[AssociationQ[result], Lookup[result, "Solver", "InvalidSolver"], "InvalidSolver"];
        resultKind = If[AssociationQ[result], Lookup[result, "ResultKind", Missing["NotAvailable"]], Missing["NotAvailable"]];
        feasibility = If[AssociationQ[result], Lookup[result, "Feasibility", Missing["NotAvailable"]], Missing["NotAvailable"]];
        message = If[AssociationQ[result], Lookup[result, "Message", Missing["NotAvailable"]], "InvalidSolverReturnShape"];
        <|
            "Method" -> method,
            "Solver" -> solver,
            "ResultKind" -> resultKind,
            "Feasibility" -> feasibility,
            "Message" -> message,
            "Decision" -> decision
        |>
    ];

SolveMFGFinalizeWithTrace[result_Association, methodUsed_, trace_List] :=
    Join[result, <|"MethodUsed" -> methodUsed, "MethodTrace" -> trace|>];

SolveMFGAbortResult[message_, methodUsed_, trace_List] :=
    <|
        "Solver" -> "SolveMFG",
        "ResultKind" -> "Failure",
        "Feasibility" -> Missing["NotApplicable"],
        "Message" -> message,
        "Method" -> "Automatic",
        "MethodUsed" -> methodUsed,
        "MethodTrace" -> trace,
        "Solution" -> Missing["NotAvailable"],
        "Convergence" -> Missing["NotApplicable"]
    |>;

SolveMFGIsCriticalQ[opts_List] :=
    Module[{alphaOpt},
        alphaOpt = "CongestionExponentFunction" /. opts;
        If[alphaOpt === "CongestionExponentFunction", alphaOpt = Automatic];
        Which[
            alphaOpt === Automatic, True,
            NumericQ[alphaOpt], alphaOpt == 1,
            AssociationQ[alphaOpt], AllTrue[Values[alphaOpt], # == 1 &],
            True, False
        ]
    ];

SolveMFGCompileInput[input_Association, opts_List:{}] :=
    Module[{potentialFunction, congestionExponentFunction, interactionFunction},
        If[SolveMFGCompiledInputQ[input], Return[input, Module]];
        potentialFunction = Replace["PotentialFunction" /. opts, "PotentialFunction" -> Automatic];
        congestionExponentFunction = Replace[
            "CongestionExponentFunction" /. opts, "CongestionExponentFunction" -> Automatic];
        interactionFunction = Replace["InteractionFunction" /. opts, "InteractionFunction" -> Automatic];
        WithHamiltonianFunctions[potentialFunction, congestionExponentFunction, interactionFunction,
            DataToEquations[input]]
    ];

SolveMFGAutomaticDispatch[input_Association, opts_List] :=
    Module[{d2e, result = Missing["NotAvailable"], decision, trace = {}, methodUsed = "Automatic",
            attempts, isCritical},
        d2e = If[SolveMFGCompiledInputQ[input], input, Quiet @ Check[SolveMFGCompileInput[input, opts], $Failed]];
        If[!AssociationQ[d2e],
            Return[SolveMFGAbortResult["DataToEquationsFailed", methodUsed, trace], Module]
        ];
        isCritical = SolveMFGIsCriticalQ[opts];
        attempts = Join[
            (* alpha=1: try CriticalCongestion first; alpha!=1: skip it *)
            If[isCritical,
                {<|
                    "Method" -> "CriticalCongestion",
                    "Run" -> Function[
                        CriticalCongestionSolver[
                            d2e,
                            Sequence @@ FilterRules[opts, Options[CriticalCongestionSolver]]
                        ]
                    ]
                |>},
                {}
            ],
            {
                <|
                    "Method" -> "Monotone",
                    "Run" -> Function[
                        MonotoneSolver[
                            d2e,
                            Sequence @@ FilterRules[opts, Options[MonotoneSolver]]
                        ]
                    ]
                |>,
                <|
                    "Method" -> "NonLinear",
                    "Run" -> Function[
                        NonLinearSolver[
                            d2e,
                            Sequence @@ FilterRules[opts, Options[NonLinearSolver]]
                        ]
                    ]
                |>
            }
        ];
        Do[
            result = attempt["Run"][];
            decision = SolveMFGClassifyOutcome[result];
            methodUsed = attempt["Method"];
            trace = Append[trace, SolveMFGBuildTraceEntry[methodUsed, result, decision]];
            Switch[decision,
                "Done",
                    Return[SolveMFGFinalizeWithTrace[result, methodUsed, trace], Module],
                "Fallback",
                    Null,
                _,
                    Return[SolveMFGAbortResult["InvalidSolverReturnShape", methodUsed, trace], Module]
            ],
            {attempt, attempts}
        ];
        If[AssociationQ[result],
            SolveMFGFinalizeWithTrace[result, methodUsed, trace],
            SolveMFGAbortResult["InvalidSolverReturnShape", methodUsed, trace]
        ]
    ];

SolveMFG[input_Association, opts:OptionsPattern[]] :=
    Module[{method, normalizedMethod, d2e},
        method = OptionValue[Method];
        normalizedMethod = ToString[method];

        Switch[normalizedMethod,
            "Monotone" | "MonotoneSolver",
                If[
                    SolveMFGCompiledInputQ[input],
                    MonotoneSolver[
                        input,
                        Sequence @@ FilterRules[{opts}, Options[MonotoneSolver]]
                    ],
                    MonotoneSolverFromData[
                        input,
                        Sequence @@ FilterRules[{opts}, Options[MonotoneSolverFromData]]
                    ]
                ]
            ,
            "CriticalCongestion" | "CriticalCongestionSolver",
                d2e = If[SolveMFGCompiledInputQ[input], input, SolveMFGCompileInput[input, {opts}]];
                CriticalCongestionSolver[
                    d2e,
                    Sequence @@ FilterRules[{opts}, Options[CriticalCongestionSolver]]
                ]
            ,
            "Automatic",
                SolveMFGAutomaticDispatch[input, {opts}]
            ,
            "NonLinear" | "NonLinearSolver",
                d2e = If[SolveMFGCompiledInputQ[input], input, SolveMFGCompileInput[input, {opts}]];
                NonLinearSolver[
                    d2e,
                    Sequence @@ FilterRules[{opts}, Options[NonLinearSolver]]
                ]
            ,
            _,
                SolveMFGUnknownMethodResult[method]
        ]
    ];

SolveMFG[_, opts:OptionsPattern[]] :=
    SolveMFGUnknownMethodResult[OptionValue[Method]];

EndPackage[];
