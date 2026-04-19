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
infeasibility or failure. When \"CongestionExponentFunction\" is known to be non-1 or \
cannot be proven to be 1, CriticalCongestion is skipped (it only handles alpha=1).
  \"CriticalCongestion\" — requires \"CongestionExponentFunction\" -> 1. Unsupported \
exponents return an immediate standardized failure without compiling raw input.
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

BuildFeasibleFlowSeed[backendState_Association, tol_: 10^-6] :=
    Module[{eqs, staticData, decoupling, massConstraints, allFlowVars, jVars, jtVars,
      massVars, aEq, bEqRaw, bEq, dims, qpVars, constraints, rawSolution, rawFlowVec,
      massAssoc, massVec, flowAssoc, spatialAssoc, transitionAssoc, flowVec, residual, nnViolation,
      comparisonData, comparableFlowVec, reason},
        eqs = Lookup[backendState, "Eqs", Missing["NotAvailable"]];
        staticData = Lookup[backendState, "StaticData", Missing["NotAvailable"]];
        decoupling = Lookup[
            eqs,
            "CriticalDecoupling",
            Lookup[Lookup[eqs, "NumericState", <||>], "CriticalDecoupling", Missing["NotAvailable"]]
        ];
        If[!AssociationQ[eqs] || !AssociationQ[staticData] || !AssociationQ[decoupling],
            Return[
                Failure["FictitiousPlaySeed", <|"Reason" -> "MissingDecouplingData"|>],
                Module
            ]
        ];

        massConstraints = Lookup[decoupling, "MassConstraints", Missing["NotAvailable"]];
        If[!AssociationQ[massConstraints],
            Return[
                Failure["FictitiousPlaySeed", <|"Reason" -> "MissingDecouplingData"|>],
                Module
            ]
        ];
        allFlowVars = Lookup[staticData, "FlowVariables", {}];
        jVars = Lookup[staticData, "JVars", {}];
        jtVars = Lookup[staticData, "JTVars", {}];
        massVars = Lookup[massConstraints, "Variables", {}];
        aEq = Lookup[massConstraints, "Aeq", Missing["NotAvailable"]];
        bEqRaw = Lookup[massConstraints, "beq", Missing["NotAvailable"]];
        If[Head[aEq] =!= SparseArray,
            aEq = Quiet @ Check[SparseArray[aEq], Missing["NotAvailable"]]
        ];
        bEq = Developer`ToPackedArray @ N @ Which[
            Head[bEqRaw] === SparseArray, Normal[bEqRaw],
            VectorQ[bEqRaw], bEqRaw,
            True, {}
        ];
        dims =
            If[Head[aEq] === SparseArray && Length[Dimensions[aEq]] === 2,
                Dimensions[aEq],
                {}
            ];
        reason =
            Which[
                !ListQ[allFlowVars] || allFlowVars =!= Join[jVars, jtVars], "InvalidStaticData",
                Complement[massVars, allFlowVars] =!= {}, "MassConstraintScopeMismatch",
                dims === {} || Length[dims] =!= 2, "InvalidMassConstraintShape",
                dims[[2]] =!= Length[massVars] || dims[[1]] =!= Length[bEq], "InvalidMassConstraintShape",
                !VectorQ[bEq, NumericQ], "InvalidMassConstraintShape",
                True, None
            ];
        If[StringQ[reason],
            Return[Failure["FictitiousPlaySeed", <|"Reason" -> reason|>], Module]
        ];

        qpVars = Table[Unique["seedVar"], {Length[massVars]}];
        constraints = Join[
            Thread[(aEq . qpVars) == bEq],
            Thread[qpVars >= 0]
        ];
        rawSolution = Quiet @ Check[
            LinearOptimization[
                Total[qpVars],
                constraints,
                qpVars
            ],
            $Failed
        ];
        If[!MatchQ[rawSolution, {_Rule ..}],
            Return[Failure["FictitiousPlaySeed", <|"Reason" -> "NoFeasibleSeed"|>], Module]
        ];
        rawFlowVec = Lookup[Association[rawSolution], qpVars, Missing["NotAvailable"]];
        If[!VectorQ[rawFlowVec, NumericQ] || Length[rawFlowVec] =!= Length[massVars],
            Return[Failure["FictitiousPlaySeed", <|"Reason" -> "NoFeasibleSeed"|>], Module]
        ];

        massAssoc = AssociationThread[massVars, N @ rawFlowVec];
        massVec = Developer`ToPackedArray @ N @ Lookup[massAssoc, massVars, Missing["NotAvailable"]];
        If[!VectorQ[massVec, NumericQ] || Length[massVec] =!= Length[massVars],
            Return[Failure["FictitiousPlaySeed", <|"Reason" -> "NoFeasibleSeed"|>], Module]
        ];
        residual = Quiet @ Check[N @ Norm[aEq . massVec - bEq, Infinity], Infinity];
        nnViolation = Quiet @ Check[N @ Max[0., -Min[massVec]], Infinity];
        If[!NumericQ[residual] || !NumericQ[nnViolation] || residual > tol || nnViolation > tol,
            Return[Failure["FictitiousPlaySeed", <|"Reason" -> "NoFeasibleSeed"|>], Module]
        ];

        flowAssoc = RecoverCriticalFlowAssociation[eqs, massAssoc, tol];
        If[MatchQ[flowAssoc, Failure[_, _Association]],
            Return[flowAssoc, Module]
        ];
        flowVec = Developer`ToPackedArray @ N @ Lookup[flowAssoc, allFlowVars, Missing["NotAvailable"]];
        If[!VectorQ[flowVec, NumericQ] || Length[flowVec] =!= Length[allFlowVars],
            Return[Failure["FictitiousPlaySeed", <|"Reason" -> "NoFeasibleSeed"|>], Module]
        ];
        spatialAssoc = KeyTake[flowAssoc, jVars];
        transitionAssoc = KeyTake[flowAssoc, jtVars];

        comparisonData = BuildSolverComparisonData[eqs, flowAssoc];
        comparableFlowVec = Lookup[comparisonData, "ComparableFlowVector", Missing["NotAvailable"]];
        If[
            !ListQ[comparableFlowVec] ||
            !VectorQ[comparableFlowVec, NumericQ] ||
            Length[comparableFlowVec] =!= Length[Lookup[eqs, "edgeList", {}]],
            Return[
                Failure["FictitiousPlaySeed", <|"Reason" -> "InvalidComparableFlowVector"|>],
                Module
            ]
        ];

        <|
            "FlowAssociation" -> flowAssoc,
            "FlowVector" -> flowVec,
            "SpatialFlowAssociation" -> spatialAssoc,
            "TransitionFlowAssociation" -> transitionAssoc,
            "NodeMassAssociation" -> Missing["NotAvailable"],
            "ComparableFlowVector" -> comparableFlowVec,
            "FeasibleQ" -> True,
            "Residuals" -> <|
                "KirchhoffResidual" -> residual,
                "NonnegativityViolation" -> nnViolation
            |>
        |>
    ];

ExtractBellmanPotentials[backendState_Association, flowState_Association, tol_: 10^-6] :=
    Module[{eqs, staticData, flowAssoc, flowVars, us, edgeList, halfPairs, switching,
      comparisonData, signedFlowVec, pairCosts, statePairs, nodePotentials, valueSystem,
      utilityAssoc, reducedCostAssoc, edgeCostAssoc, bellmanResidual, allCosts,
      utilityVals, reducedVals, reason},
        eqs = Lookup[backendState, "Eqs", Missing["NotAvailable"]];
        staticData = Lookup[backendState, "StaticData", Missing["NotAvailable"]];
        flowAssoc = Lookup[flowState, "FlowAssociation", Missing["NotAvailable"]];
        flowVars = Lookup[staticData, "FlowVariables", {}];
        us = Lookup[eqs, "us", {}];
        edgeList = Lookup[eqs, "edgeList", {}];
        halfPairs = List @@@ edgeList;
        switching = Lookup[eqs, "SwitchingCosts", <||>];
        reason =
            Which[
                !AssociationQ[eqs] || !AssociationQ[staticData], "MissingBackendState",
                !AssociationQ[flowAssoc], "MissingFlowState",
                !ListQ[flowVars] || !And @@ (KeyExistsQ[flowAssoc, #] & /@ flowVars), "IncompleteFlowState",
                !ListQ[us], "InvalidUtilityState",
                True, None
            ];
        If[StringQ[reason],
            Return[Failure["FictitiousPlayPotentials", <|"Reason" -> reason|>], Module]
        ];

        comparisonData = BuildSolverComparisonData[eqs, flowAssoc];
        signedFlowVec = Lookup[comparisonData, "ComparableFlowVector", Missing["NotAvailable"]];
        If[
            !ListQ[signedFlowVec] ||
            !VectorQ[signedFlowVec, NumericQ] ||
            Length[signedFlowVec] =!= Length[edgeList],
            Return[Failure["FictitiousPlayPotentials", <|"Reason" -> "MissingComparableFlowVector"|>], Module]
        ];

        pairCosts =
            If[halfPairs === {},
                <||>,
                BuildMonotonePairCostAssociation[halfPairs, edgeList, signedFlowVec]
            ];
        edgeCostAssoc = Join[
            Association @ Map[
                Function[{var},
                    var -> N @ LookupAssociationValue[pairCosts, List @@ var, 0.]
                ],
                Lookup[eqs, "js", {}]
            ],
            Association @ Map[
                Function[{var},
                    var -> N @ (
                        LookupAssociationValue[switching, List @@ var, 0.] +
                        LookupAssociationValue[pairCosts, List @@ var[[{2, 3}]], 0.]
                    )
                ],
                Lookup[eqs, "jts", {}]
            ]
        ];
        allCosts = Values[edgeCostAssoc];
        If[!VectorQ[allCosts, NumericQ],
            Return[Failure["FictitiousPlayPotentials", <|"Reason" -> "NonNumericEdgeCost"|>], Module]
        ];
        If[allCosts =!= {} && Min[allCosts] < -tol,
            Return[Failure["FictitiousPlayPotentials", <|"Reason" -> "NegativeEdgeCost"|>], Module]
        ];

        valueSystem = BuildMonotoneValueSystem[eqs];
        utilityAssoc = Lookup[valueSystem, "StateValueAssociation", Missing["NotAvailable"]][pairCosts];
        If[!AssociationQ[utilityAssoc],
            Return[Failure["FictitiousPlayPotentials", <|"Reason" -> "BellmanSolveFailed"|>], Module]
        ];
        utilityVals = Lookup[utilityAssoc, us, Missing["NotAvailable"]];
        If[!VectorQ[utilityVals, NumericQ] || Length[utilityVals] =!= Length[us],
            Return[Failure["FictitiousPlayPotentials", <|"Reason" -> "BellmanSolveFailed"|>], Module]
        ];
        utilityAssoc = AssociationThread[us, N @ utilityVals];
        statePairs = DeleteDuplicates @ Cases[us, u[a_, b_] :> {a, b}];
        nodePotentials = AssociationThread[statePairs, Lookup[utilityAssoc, u @@@ statePairs]];

        reducedCostAssoc = Join[
            Association @ Map[
                Function[{var},
                    With[{pair = List @@ var},
                        var -> N @ (
                            LookupAssociationValue[pairCosts, pair, 0.] +
                            LookupAssociationValue[utilityAssoc, u @@ pair, 0.] -
                            LookupAssociationValue[utilityAssoc, u @@ Reverse[pair], 0.]
                        )
                    ]
                ],
                Lookup[eqs, "js", {}]
            ],
            Association @ Map[
                Function[{var},
                    var -> N @ (
                        LookupAssociationValue[switching, List @@ var, 0.] +
                        LookupAssociationValue[pairCosts, List @@ var[[{2, 3}]], 0.] +
                        LookupAssociationValue[utilityAssoc, u @@ (List @@ var[[{2, 3}]]), 0.] -
                        LookupAssociationValue[utilityAssoc, u @@ (List @@ var[[{1, 2}]]), 0.]
                    )
                ],
                Lookup[eqs, "jts", {}]
            ]
        ];
        reducedVals = Lookup[reducedCostAssoc, flowVars, Missing["NotAvailable"]];
        If[!VectorQ[reducedVals, NumericQ] || Length[reducedVals] =!= Length[flowVars],
            Return[Failure["FictitiousPlayPotentials", <|"Reason" -> "InvalidReducedCosts"|>], Module]
        ];
        reducedCostAssoc = AssociationThread[flowVars, N @ reducedVals];
        bellmanResidual =
            With[{transitionReducedVals = Lookup[reducedCostAssoc, Lookup[eqs, "jts", {}], {}]},
            If[transitionReducedVals === {},
                0.,
                N @ Max[0., -Min[transitionReducedVals]]
            ]];
        If[!NumericQ[bellmanResidual] || bellmanResidual > tol,
            Return[Failure["FictitiousPlayPotentials", <|"Reason" -> "BellmanResidualViolation"|>], Module]
        ];

        <|
            "NodePotentialAssociation" -> nodePotentials,
            "UtilityAssociation" -> utilityAssoc,
            "ReducedCostAssociation" -> reducedCostAssoc,
            "EdgeCostAssociation" -> edgeCostAssoc,
            "BellmanResidual" -> bellmanResidual,
            "CostNonnegativeQ" -> True
        |>
    ];

Options[BuildSoftPolicyAndPropagate] = {
    "Temperature" -> 0.1,
    "Damping" -> 0.5
};

BuildSoftPolicyAndPropagate[
    backendState_Association,
    flowState_Association,
    potentialState_Association,
    opts : OptionsPattern[]
] :=
    Module[{eqs, staticData, decoupling, massConstraints, allFlowVars, jVars, jtVars,
      oldFlowAssoc, reducedCostAssoc, auxTriples, us, statePairs, outgoingByState,
      eta, alpha, topoData, stateGraph, topoOrder, policyByNode, flatPolicyWeights = <||>,
      nodeMass, entryRules, entryPairs, state, outgoingTriples, outgoingVars, costs,
      minCost, shiftedCosts, weights, totalWeight, probs, mass, target, flow, jtsAssoc,
      jsAssoc, partialPushedAssoc, pushedFlowAssoc, oldFlowVec, pushedFlowVec, newFlowVec,
      newFlowAssoc, spatialAssoc, transitionAssoc, comparableFlowVec, comparisonData,
      massVars, aEq, bEqRaw, bEq, massVec, residual, nnViolation, newNodeMassAssoc, reason},
        eqs = Lookup[backendState, "Eqs", Missing["NotAvailable"]];
        staticData = Lookup[backendState, "StaticData", Missing["NotAvailable"]];
        oldFlowAssoc = Lookup[flowState, "FlowAssociation", Missing["NotAvailable"]];
        reducedCostAssoc = Lookup[potentialState, "ReducedCostAssociation", Missing["NotAvailable"]];
        eta = N @ OptionValue["Temperature"];
        alpha = N @ OptionValue["Damping"];
        decoupling = Lookup[
            eqs,
            "CriticalDecoupling",
            Lookup[Lookup[eqs, "NumericState", <||>], "CriticalDecoupling", Missing["NotAvailable"]]
        ];
        massConstraints = Lookup[decoupling, "MassConstraints", Missing["NotAvailable"]];
        allFlowVars = Lookup[staticData, "FlowVariables", {}];
        jVars = Lookup[staticData, "JVars", Lookup[eqs, "js", {}]];
        jtVars = Lookup[staticData, "JTVars", Lookup[eqs, "jts", {}]];
        auxTriples = Lookup[eqs, "auxTriples", {}];
        us = Lookup[eqs, "us", {}];
        reason =
            Which[
                !AssociationQ[eqs] || !AssociationQ[staticData], "MissingBackendState",
                !AssociationQ[oldFlowAssoc], "MissingFlowState",
                !AssociationQ[reducedCostAssoc], "MissingPotentialState",
                !And @@ (KeyExistsQ[oldFlowAssoc, #] & /@ allFlowVars), "IncompleteFlowState",
                !And @@ (KeyExistsQ[reducedCostAssoc, #] & /@ allFlowVars), "IncompletePotentialState",
                !NumericQ[eta] || eta <= 0, "InvalidPolicyHyperparameters",
                !NumericQ[alpha] || alpha < 0 || alpha > 1, "InvalidPolicyHyperparameters",
                !AssociationQ[decoupling] || !AssociationQ[massConstraints], "MissingDecouplingData",
                True, None
            ];
        If[StringQ[reason],
            Return[Failure["FictitiousPlayPolicy", <|"Reason" -> reason|>], Module]
        ];

        topoData = Lookup[staticData, "TopologicalData", Missing["NotAvailable"]];
        If[AssociationQ[topoData] && TrueQ[Lookup[topoData, "IsDAG", False]] === False,
            Return[Failure["FictitiousPlayPolicy", <|"Reason" -> "NonDAGSupport"|>], Module]
        ];
        statePairs = DeleteDuplicates @ Cases[us, u[a_, b_] :> {a, b}];
        stateGraph = Graph[
            statePairs,
            (DirectedEdge[#[[{1, 2}]], #[[{2, 3}]]] & /@ auxTriples),
            DirectedEdges -> True
        ];
        topoOrder = Quiet @ Check[TopologicalSort[stateGraph], {}];
        If[!ListQ[topoOrder] || topoOrder === {} || !AcyclicGraphQ[stateGraph],
            Return[Failure["FictitiousPlayPolicy", <|"Reason" -> "NonDAGSupport"|>], Module]
        ];

        outgoingByState = GroupBy[auxTriples, #[[{1, 2}]] &];
        policyByNode =
            Association @ KeyValueMap[
                Function[{sourceState, triples},
                    outgoingVars = j @@@ triples;
                    costs = Lookup[reducedCostAssoc, outgoingVars, Missing["NotAvailable"]];
                    If[VectorQ[costs, NumericQ] && outgoingVars =!= {},
                        minCost = Min[costs];
                        shiftedCosts = costs - minCost;
                        weights = Quiet[Exp[-shiftedCosts/eta], General::munfl];
                        totalWeight = Total[weights];
                        probs =
                            If[!NumericQ[totalWeight] || totalWeight <= 0,
                                ConstantArray[1./Length[outgoingVars], Length[outgoingVars]],
                                N @ (weights/totalWeight)
                            ];
                        AssociateTo[flatPolicyWeights, AssociationThread[outgoingVars, probs]];
                        sourceState -> AssociationThread[outgoingVars, probs],
                        Nothing
                    ]
                ],
                outgoingByState
            ];

        nodeMass = AssociationThread[statePairs, ConstantArray[0., Length[statePairs]]];
        entryRules = Association @ Flatten[ToRules /@ Lookup[eqs, "EqEntryIn", {}]];
        entryPairs = Select[Lookup[eqs, "js", {}], KeyExistsQ[entryRules, #] &];
        Do[
            AssociateTo[nodeMass, List @@ pairVar -> N @ entryRules[pairVar]],
            {pairVar, entryPairs}
        ];

        jtsAssoc = AssociationThread[jtVars, ConstantArray[0., Length[jtVars]]];
        Do[
            state = topoOrder[[k]];
            mass = N @ LookupAssociationValue[nodeMass, state, 0.];
            outgoingTriples = LookupAssociationValue[outgoingByState, state, {}];
            If[outgoingTriples === {} || !NumericQ[mass] || mass == 0.,
                Continue[]
            ];
            Do[
                flow = N @ (mass * LookupAssociationValue[flatPolicyWeights, j @@ triple, 0.]);
                jtsAssoc[j @@ triple] = flow;
                target = triple[[{2, 3}]];
                nodeMass[target] = N @ (LookupAssociationValue[nodeMass, target, 0.] + flow),
                {triple, outgoingTriples}
            ],
            {k, Length[topoOrder]}
        ];

        jsAssoc = AssociationThread[
            jVars,
            N @ (Lookup[nodeMass, List @@@ jVars, 0.])
        ];
        partialPushedAssoc = Join[jsAssoc, jtsAssoc];
        pushedFlowAssoc = RecoverCriticalFlowAssociation[eqs, partialPushedAssoc, 10^-6];
        If[MatchQ[pushedFlowAssoc, Failure[_, _Association]],
            Return[Failure["FictitiousPlayPolicy", <|"Reason" -> "PropagationInfeasible"|>], Module]
        ];

        oldFlowVec = Lookup[oldFlowAssoc, allFlowVars, Missing["NotAvailable"]];
        pushedFlowVec = Lookup[pushedFlowAssoc, allFlowVars, Missing["NotAvailable"]];
        If[
            !VectorQ[oldFlowVec, NumericQ] ||
            !VectorQ[pushedFlowVec, NumericQ] ||
            Length[oldFlowVec] =!= Length[allFlowVars] ||
            Length[pushedFlowVec] =!= Length[allFlowVars],
            Return[Failure["FictitiousPlayPolicy", <|"Reason" -> "PropagationInfeasible"|>], Module]
        ];
        newFlowVec = N @ ((1. - alpha) oldFlowVec + alpha pushedFlowVec);
        newFlowAssoc = AssociationThread[allFlowVars, newFlowVec];

        massVars = Lookup[massConstraints, "Variables", {}];
        aEq = Lookup[massConstraints, "Aeq", Missing["NotAvailable"]];
        bEqRaw = Lookup[massConstraints, "beq", Missing["NotAvailable"]];
        If[Head[aEq] =!= SparseArray,
            aEq = Quiet @ Check[SparseArray[aEq], Missing["NotAvailable"]]
        ];
        bEq = Developer`ToPackedArray @ N @ Which[
            Head[bEqRaw] === SparseArray, Normal[bEqRaw],
            VectorQ[bEqRaw], bEqRaw,
            True, {}
        ];
        massVec = Lookup[newFlowAssoc, massVars, Missing["NotAvailable"]];
        residual = Quiet @ Check[N @ Norm[aEq . massVec - bEq, Infinity], Infinity];
        nnViolation = Quiet @ Check[N @ Max[0., -Min[newFlowVec]], Infinity];
        If[
            !VectorQ[massVec, NumericQ] ||
            !NumericQ[residual] ||
            !NumericQ[nnViolation] ||
            residual > 10^-6 ||
            nnViolation > 10^-6,
            Return[Failure["FictitiousPlayPolicy", <|"Reason" -> "PropagationInfeasible"|>], Module]
        ];

        spatialAssoc = KeyTake[newFlowAssoc, jVars];
        transitionAssoc = KeyTake[newFlowAssoc, jtVars];
        comparisonData = BuildSolverComparisonData[eqs, newFlowAssoc];
        comparableFlowVec = Lookup[comparisonData, "ComparableFlowVector", Missing["NotAvailable"]];
        newNodeMassAssoc = AssociationThread[List @@@ jVars, Lookup[newFlowAssoc, jVars, 0.]];
        If[
            !ListQ[comparableFlowVec] ||
            !VectorQ[comparableFlowVec, NumericQ] ||
            Length[comparableFlowVec] =!= Length[Lookup[eqs, "edgeList", {}]],
            Return[Failure["FictitiousPlayPolicy", <|"Reason" -> "PropagationInfeasible"|>], Module]
        ];

        <|
            "PolicyState" -> <|
                "OutgoingPolicyByNode" -> policyByNode,
                "OutgoingArcWeights" -> flatPolicyWeights,
                "Temperature" -> eta,
                "Damping" -> alpha,
                "DirectedSupportGraph" -> stateGraph,
                "TopologicalOrder" -> topoOrder,
                "PropagationMode" -> "DAGForwardPass"
            |>,
            "FlowState" -> <|
                "FlowAssociation" -> newFlowAssoc,
                "FlowVector" -> Developer`ToPackedArray @ N @ newFlowVec,
                "SpatialFlowAssociation" -> spatialAssoc,
                "TransitionFlowAssociation" -> transitionAssoc,
                "NodeMassAssociation" -> newNodeMassAssoc,
                "ComparableFlowVector" -> comparableFlowVec,
                "FeasibleQ" -> True,
                "Residuals" -> <|
                    "KirchhoffResidual" -> residual,
                    "NonnegativityViolation" -> nnViolation
                |>
            |>
        |>
    ];

Options[ClassifyAndCheckStability] = {
    "FlowThresholdOn" -> 10^-4,
    "FlowThresholdOff" -> 10^-5,
    "ReducedCostEq" -> 10^-4,
    "ReducedCostGap" -> 10^-3,
    "StableIterationsLimit" -> 3
};

ClassifyAndCheckStability[
    backendState_Association,
    flowState_Association,
    potentialState_Association,
    history_Association : <||>,
    opts : OptionsPattern[]
] :=
    Module[{staticData, flowAssoc, reducedCostAssoc, flowVars, epsOn, epsOff, deltaEq, deltaGap,
      stableLimit, activeVars, inactiveVars, ambiguousVars, currentSignature, previousSignature,
      previousStableCount, currentStableCount, stableQ, supportChangedQ, allCoveredQ},
        staticData = Lookup[backendState, "StaticData", Missing["NotAvailable"]];
        flowAssoc = Lookup[flowState, "FlowAssociation", Missing["NotAvailable"]];
        reducedCostAssoc = Lookup[potentialState, "ReducedCostAssociation", Missing["NotAvailable"]];
        flowVars = Lookup[staticData, "FlowVariables", {}];
        epsOn = N @ OptionValue["FlowThresholdOn"];
        epsOff = N @ OptionValue["FlowThresholdOff"];
        deltaEq = N @ OptionValue["ReducedCostEq"];
        deltaGap = N @ OptionValue["ReducedCostGap"];
        stableLimit = OptionValue["StableIterationsLimit"];
        If[
            !AssociationQ[staticData] ||
            !AssociationQ[flowAssoc] ||
            !AssociationQ[reducedCostAssoc] ||
            !ListQ[flowVars] ||
            !And @@ (KeyExistsQ[flowAssoc, #] & /@ flowVars) ||
            !And @@ (KeyExistsQ[reducedCostAssoc, #] & /@ flowVars) ||
            !NumericQ[epsOn] || !NumericQ[epsOff] || !NumericQ[deltaEq] || !NumericQ[deltaGap] ||
            !IntegerQ[stableLimit] || stableLimit < 1 ||
            epsOff > epsOn,
            Return[Failure["FictitiousPlayClassification", <|"Reason" -> "InvalidClassificationInput"|>], Module]
        ];

        activeVars = Select[
            flowVars,
            Lookup[flowAssoc, #, -Infinity] >= epsOn &&
            Lookup[reducedCostAssoc, #, Infinity] <= deltaEq &
        ];
        inactiveVars = Select[
            Complement[flowVars, activeVars],
            Lookup[flowAssoc, #, Infinity] <= epsOff &&
            Lookup[reducedCostAssoc, #, -Infinity] >= deltaGap &
        ];
        ambiguousVars = Complement[flowVars, Join[activeVars, inactiveVars]];
        allCoveredQ =
            Sort[Join[activeVars, inactiveVars, ambiguousVars]] === Sort[flowVars] &&
            Intersection[activeVars, inactiveVars] === {} &&
            Intersection[activeVars, ambiguousVars] === {} &&
            Intersection[inactiveVars, ambiguousVars] === {};
        If[!allCoveredQ,
            Return[Failure["FictitiousPlayClassification", <|"Reason" -> "InvalidClassificationPartition"|>], Module]
        ];

        currentSignature = <|
            "Active" -> SortBy[activeVars, ToString[#, InputForm] &],
            "Inactive" -> SortBy[inactiveVars, ToString[#, InputForm] &]
        |>;
        previousSignature = Lookup[history, "LastSupportSignature", None];
        previousStableCount = Lookup[history, "StableIterations", 0];
        supportChangedQ = !(AssociationQ[previousSignature] && previousSignature === currentSignature);
        currentStableCount =
            If[supportChangedQ,
                0,
                previousStableCount + 1
            ];
        stableQ = currentStableCount >= stableLimit;

        <|
            "ClassificationState" -> <|
                "ActiveVariables" -> activeVars,
                "InactiveVariables" -> inactiveVars,
                "AmbiguousVariables" -> ambiguousVars,
                "FlowThresholds" -> <|"On" -> epsOn, "Off" -> epsOff|>,
                "ReducedCostThresholds" -> <|"Eq" -> deltaEq, "Gap" -> deltaGap|>,
                "StableQ" -> stableQ,
                "StableIterations" -> currentStableCount,
                "SupportSignature" -> currentSignature,
                "SupportChangedQ" -> supportChangedQ
            |>,
            "OracleState" -> <|
                "OracleReadyQ" -> stableQ,
                "PrunedEqualities" -> If[stableQ, activeVars, Missing["NotAvailable"]],
                "PrunedZeroFlows" -> If[stableQ, inactiveVars, Missing["NotAvailable"]],
                "ResidualDisjunctions" -> If[stableQ, ambiguousVars, Missing["NotAvailable"]],
                "HandoffReason" -> If[stableQ, "StableSupport", "NotReady"]
            |>
        |>
    ];

(* INTERNAL ONLY — SolveCriticalFictitiousPlayBackend

   Fictitious Play backend for critical congestion support discovery.

   *** KNOWN VULNERABILITIES (do NOT expose publicly until fixed) ***
   1. Cost-scale fragility: softmax temperature assumes O(1) costs; fails on wildly different scales
   2. Symmetric graph oscillation: tied best-response edges prevent convergence, dump to symbolic solver
   3. Float-to-exact precision chasm: tolerance-based zero classification (10^-5) can over-constrain symbolic EE
   4. DAG cycle brittleness: no recovery from floating-point-induced cycles; instant abort

   See project memory "Phase 5/6 Technical Debt" for detailed mitigations and integration checklist.

   v1 constraints: DAG-only, cost scales O(1)–O(100), avoid highly symmetric topologies.
   Status: internal use only; not wired into public SolveMFG API. *)

Options[SolveCriticalFictitiousPlayBackend] = {
    "MaxIterations" -> 20,
    "Temperature" -> 0.1,
    "Damping" -> 0.5,
    "FlowThresholdOn" -> 10^-4,
    "FlowThresholdOff" -> 10^-5,
    "ReducedCostEq" -> 10^-4,
    "ReducedCostGap" -> 10^-3,
    "StableIterationsLimit" -> 3
};

SolveCriticalFictitiousPlayBackend[
    backendState_Association,
    opts : OptionsPattern[]
] :=
    Module[{
        seedFlowState, maxIter, allStates, finalState, iterationLog,
        stopReason, resultKind, message
    },
        (* Initialize seed flow *)
        seedFlowState = BuildFeasibleFlowSeed[backendState];

        (* Check seed feasibility *)
        If[!Lookup[seedFlowState, "FeasibleQ", False],
            Return[<|
                "Solver" -> "CriticalFictitiousPlay",
                "ResultKind" -> "Failure",
                "Feasibility" -> Missing["NotAvailable"],
                "Message" -> "SeedInfeasible",
                "Solution" -> Missing["NotAvailable"],
                "FlowState" -> seedFlowState,
                "PotentialState" -> Missing["NotAvailable"],
                "PolicyState" -> Missing["NotAvailable"],
                "ClassificationState" -> Missing["NotAvailable"],
                "OracleState" -> <|"OracleReadyQ" -> False|>,
                "History" -> <|
                    "StopReason" -> "SeedInfeasibility",
                    "IterationCount" -> 0,
                    "LastSupportSignature" -> Missing["NotAvailable"],
                    "IterationLog" -> {}
                |>,
                "Diagnostics" -> <|"SeedFeasible" -> False|>
            |>]
        ];

        (* Extract options *)
        maxIter = OptionValue["MaxIterations"];

        (* NestWhileList step function over compound state *)
        allStates = NestWhileList[
            Function[state,
                Module[{potState, policyResult, classResult},
                    potState = ExtractBellmanPotentials[backendState, state["FlowState"]];
                    policyResult = BuildSoftPolicyAndPropagate[
                        backendState,
                        state["FlowState"],
                        potState,
                        Sequence @@ FilterRules[{opts}, Options[BuildSoftPolicyAndPropagate]]
                    ];
                    classResult = ClassifyAndCheckStability[
                        backendState,
                        policyResult["FlowState"],
                        potState,
                        state["ThinHistory"],
                        Sequence @@ FilterRules[{opts}, Options[ClassifyAndCheckStability]]
                    ];
                    <|
                        "Iteration" -> state["Iteration"] + 1,
                        "FlowState" -> policyResult["FlowState"],
                        "PotentialState" -> potState,
                        "PolicyState" -> policyResult["PolicyState"],
                        "ClassificationState" -> classResult["ClassificationState"],
                        "OracleState" -> classResult["OracleState"],
                        "ThinHistory" -> <|
                            "LastSupportSignature" -> classResult["ClassificationState"]["SupportSignature"],
                            "StableIterations" -> classResult["ClassificationState"]["StableIterations"]
                        |>
                    |>
                ]
            ],
            <|
                "Iteration" -> 0,
                "FlowState" -> seedFlowState,
                "PotentialState" -> Missing["NotAvailable"],
                "PolicyState" -> Missing["NotAvailable"],
                "ClassificationState" -> Missing["NotAvailable"],
                "OracleState" -> <|"OracleReadyQ" -> False|>,
                "ThinHistory" -> <||>
            |>,
            Not[#["OracleState"]["OracleReadyQ"]] &,
            1,
            maxIter
        ];

        finalState = Last[allStates];

        (* Build iteration log from allStates[[2;;]] *)
        iterationLog = If[Length[allStates] > 1,
            Table[
                With[{s = allStates[[k]]},
                    <|
                        "Iteration" -> s["Iteration"],
                        "BellmanResidual" -> Lookup[s["PotentialState"], "BellmanResidual", Missing["NotAvailable"]],
                        "SupportSignature" -> Lookup[s["ClassificationState"], "SupportSignature", Missing["NotAvailable"]],
                        "StableIterations" -> Lookup[s["ClassificationState"], "StableIterations", Missing["NotAvailable"]],
                        "OracleReadyQ" -> s["OracleState"]["OracleReadyQ"]
                    |>
                ],
                {k, 2, Length[allStates]}
            ],
            {}
        ];

        (* Determine stop reason and result kind *)
        If[finalState["OracleState"]["OracleReadyQ"],
            stopReason = "StableSupport";
            resultKind = "Success";
            message = None,
            stopReason = "MaxIterations";
            resultKind = "NonConverged";
            message = "StableSupportNotFound"
        ];

        (* Return standardized backend result *)
        <|
            "Solver" -> "CriticalFictitiousPlay",
            "ResultKind" -> resultKind,
            "Feasibility" -> If[resultKind === "Success", "Feasible", Missing["NotAvailable"]],
            "Message" -> message,
            "Solution" -> If[resultKind === "Success",
                finalState["FlowState"]["FlowAssociation"],
                Missing["NotAvailable"]
            ],
            "FlowState" -> finalState["FlowState"],
            "PotentialState" -> finalState["PotentialState"],
            "PolicyState" -> finalState["PolicyState"],
            "ClassificationState" -> finalState["ClassificationState"],
            "OracleState" -> finalState["OracleState"],
            "History" -> <|
                "StopReason" -> stopReason,
                "IterationCount" -> finalState["Iteration"],
                "LastSupportSignature" -> finalState["ClassificationState"]["SupportSignature"],
                "IterationLog" -> iterationLog
            |>,
            "Diagnostics" -> <|
                "SeedFeasible" -> True,
                "OptionsUsed" -> {opts}
            |>
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
Get["MFGraphs`Solvers`"];
Get["MFGraphs`NonLinearSolver`"];
Get["MFGraphs`Monotone`"];
Get["MFGraphs`TimeDependentSolver`"];
Get["MFGraphs`Graphics`"];

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

SolveMFGHeldOptionValue[heldOpts_HoldComplete, key_] :=
    Module[{matches},
        matches = Cases[
            heldOpts,
            HoldPattern[Rule[key, value_]] :> HoldComplete[value],
            Infinity
        ];
        If[matches === {}, HoldComplete[Automatic], Last[matches]]
    ];

SolveMFGUnsupportedCriticalAlphaResult[] :=
    <|
        "Solver" -> "CriticalCongestion",
        "ResultKind" -> "Failure",
        "Feasibility" -> Missing["NotAvailable"],
        "Message" -> "Method -> 'CriticalCongestion' only supports CongestionExponentFunction -> 1.",
        "Solution" -> Missing["NotAvailable"],
        "Convergence" -> Missing["NotApplicable"]
    |>;

SolveMFGIsCriticalQ[heldOpts_HoldComplete] :=
    Module[{heldAlphaOpt, alphaOpt},
        heldAlphaOpt = SolveMFGHeldOptionValue[heldOpts, "CongestionExponentFunction"];
        If[
            MatchQ[
                heldAlphaOpt,
                HoldComplete[
                    HoldPattern[1 &] | HoldPattern[Function[_, 1]]
                ]
            ],
            Return[True, Module]
        ];
        alphaOpt = ReleaseHold[heldAlphaOpt];
        Which[
            alphaOpt === Automatic, True,
            NumericQ[alphaOpt], alphaOpt == 1,
            AssociationQ[alphaOpt], AllTrue[Values[alphaOpt], # == 1 &],
            True, False
        ]
    ];

SolveMFGCompileInput[input_Association, opts_List:{}] :=
    Module[{potentialFunction, congestionExponentFunction, interactionFunction},
        If[SolveMFGCompiledInputQ[input],
            Return[input, Module]
        ];
        potentialFunction = Replace["PotentialFunction" /. opts, "PotentialFunction" -> Automatic];
        congestionExponentFunction = Replace[
            "CongestionExponentFunction" /. opts,
            "CongestionExponentFunction" -> Automatic
        ];
        interactionFunction = Replace["InteractionFunction" /. opts, "InteractionFunction" -> Automatic];
        WithHamiltonianFunctions[
            potentialFunction,
            congestionExponentFunction,
            interactionFunction,
            DataToEquations[input]
        ]
    ];

SolveMFGAutomaticDispatch[input_Association, opts_List] :=
    Module[{d2e, result = Missing["NotAvailable"], decision, trace = {}, methodUsed = "Automatic",
            attempts, isCritical, heldOpts},
        heldOpts = HoldComplete[{opts}];
        d2e = If[SolveMFGCompiledInputQ[input], input, Quiet @ Check[SolveMFGCompileInput[input, opts], $Failed]];
        If[!AssociationQ[d2e],
            Return[SolveMFGAbortResult["DataToEquationsFailed", methodUsed, trace], Module]
        ];
        isCritical = SolveMFGIsCriticalQ[heldOpts];
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
    Module[{method, normalizedMethod, d2e, heldOpts, isCritical},
        heldOpts = HoldComplete[{opts}];
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
                If[
                    SolveMFGCompiledInputQ[input],
                    d2e = input;
                    CriticalCongestionSolver[
                        d2e,
                        Sequence @@ FilterRules[{opts}, Options[CriticalCongestionSolver]]
                    ],
                    isCritical = SolveMFGIsCriticalQ[heldOpts];
                    If[!TrueQ[isCritical],
                        SolveMFGUnsupportedCriticalAlphaResult[],
                        d2e = SolveMFGCompileInput[input, {opts}];
                        CriticalCongestionSolver[
                            d2e,
                            Sequence @@ FilterRules[{opts}, Options[CriticalCongestionSolver]]
                        ]
                    ]
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
