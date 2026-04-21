(* Wolfram Language Package *)
(*
   MFGraphs: A Wolfram Language package for Mean Field Games on Networks.
   
   Main Components:
   - DataToEquations: Symbolic graph-to-equation converter.
   - Solvers: Critical-congestion solver suite (extracted Phase 3).
   - DNFReduce: Symbolic logical reduction engine.
   - Graphics: Public visualization helpers.
   - Scenario: Typed scenario kernel (makeScenario, validateScenario, etc.).
*)

(* Created by the Wolfram Workbench May 5, 2020 *)
(* To distribute, use ideas from https://community.wolfram.com/groups/-/m/t/214901 *)

BeginPackage["MFGraphs`"];

(* --- Public API declarations (Usage strings) --- *)

SolveMFG::usage =
"SolveMFG[input, opts] solves the critical congestion MFG and returns its \
standardized result association.

input can be either raw data (association accepted by DataToEquations) or an already \
compiled equation association.";

DataToEquations::usage =
"DataToEquations[input] converts raw network data into a compiled equation association.";

CriticalCongestionSolver::usage =
"CriticalCongestionSolver[d2e, opts] solves the critical congestion MFG system.";

MFGPreprocessing::usage =
"MFGPreprocessing[eqs] performs symbolic preprocessing and utility reduction on the MFG system.";

MFGSystemSolver::usage =
"MFGSystemSolver[eqs][approx] returns a symbolic solution for the MFG system given flow approximations.";

IsCriticalSolution::usage =
"IsCriticalSolution[eqs] validates whether a solution satisfies the critical congestion MFG equations.";

IsFeasible::usage =
"IsFeasible[result] returns True if the solver result indicates a feasible solution \
(all flows non-negative and constraints satisfied within tolerance).";

CheckFlowFeasibility::usage =
"CheckFlowFeasibility[assoc] returns \"Feasible\" if all flow variables in assoc \
are non-negative, and \"Infeasible\" otherwise.";

NumberVectorQ::usage =
"NumberVectorQ[j] returns True if the vector j is numeric.";

BuildSolverComparisonData::usage =
"BuildSolverComparisonData[Eqs, solution] returns an association with comparison \
metrics like Kirchhoff residual and boundary mass balance.";

BuildMonotoneStateData::usage =
"BuildMonotoneStateData[d2e] returns shared linear state for monotone-like solvers.";

BuildMonotoneValueSystem::usage =
"BuildMonotoneValueSystem[d2e] returns a function to solve for node potentials.";

BuildReducedKirchhoffCoordinates::usage =
"BuildReducedKirchhoffCoordinates[d2e, basePoint] returns reduced affine coordinates \
on the Kirchhoff manifold.";

BuildMonotonePairCostAssociation::usage =
"BuildMonotonePairCostAssociation[halfPairs, edgeList, q] returns signed edge costs.";

EncodeFlowAssociation::usage =
"EncodeFlowAssociation[ns, assoc] returns a packed numeric vector from a flow association.";

DecodeFlowVector::usage =
"DecodeFlowVector[ns, vec] returns a flow association from a packed numeric vector.";

ComputeSignedEdgeFlowsFast::usage =
"ComputeSignedEdgeFlowsFast[ns, vec] returns a vector of signed edge flows.";

ComputeKirchhoffResidualFast::usage =
"ComputeKirchhoffResidualFast[ns, vec] returns the maximum Kirchhoff residual.";

BuildFeasibleFlowSeed::usage =
"BuildFeasibleFlowSeed[backendState] returns an initial feasible flow for FP.";

ExtractBellmanPotentials::usage =
"ExtractBellmanPotentials[backendState, flowState] extracts node potentials.";

BuildSoftPolicyAndPropagate::usage =
"BuildSoftPolicyAndPropagate[backendState, flowState, potentialState] propagates policy.";

ClassifyAndCheckStability::usage =
"ClassifyAndCheckStability[backendState, flowState, potentialState, history] checks stability.";

SolveCriticalFictitiousPlayBackend::usage =
"SolveCriticalFictitiousPlayBackend[backendState] runs the FP iterative solver.";

SolveCriticalJFirstBackend::usage =
"SolveCriticalJFirstBackend[Eqs, opts] runs the j-first numeric strategy.";

SolveCriticalJFirstUtilities::usage =
"SolveCriticalJFirstUtilities[eqs, flowAssoc, uVars, tol] recovers node potentials from flows.";

BuildCriticalQuadraticObjective::usage =
"BuildCriticalQuadraticObjective[d2e] returns a quadratic approximation of the MFG objective.";

UseQuadraticCriticalBackendQ::usage =
"UseQuadraticCriticalBackendQ[d2e] returns True if the case is eligible for a quadratic shortcut.";

MonotoneVariableFieldValue::usage =
"MonotoneVariableFieldValue[var, values, switching] returns the field value for a variable.";

LookupAssociationValue::usage =
"LookupAssociationValue[assoc, key, default] is a robust Lookup helper.";

SelectFlowAssociation::usage =
"SelectFlowAssociation[assoc] returns a sub-association with only flow (j) keys.";

BuildBoundaryMassData::usage =
"BuildBoundaryMassData[Eqs, flowAssoc] returns mass balance metrics.";

BuildUtilityReductionResidualData::usage =
"BuildUtilityReductionResidualData[Eqs, solution] returns utility reduction metrics.";

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

MFGParallelMap::usage =
"MFGParallelMap[f, list] applies f to each element of list, using ParallelMap
when Length[list] >= $MFGraphsParallelThreshold (and kernels are already running)
or Length[list] >= $MFGraphsParallelLaunchThreshold (when kernels must be launched).";

$MFGraphsParallelLaunchThreshold::usage =
"$MFGraphsParallelLaunchThreshold is the minimum list length required to justify
launching parallel subkernels. Only applies when $KernelCount === 0. Default is 50.";

(* Load submodules in dependency order *)
Get["MFGraphs`Examples`ExamplesData`"];
Get["MFGraphs`DNFReduce`"];
Get["MFGraphs`DataToEquations`"];
Get["MFGraphs`Solvers`"];
Get["MFGraphs`Graphics`"];
Get["MFGraphs`Scenario`"];

Options[SolveMFG] = DeleteDuplicatesBy[
    Join[
        {Method -> "Automatic"},
        Options[CriticalCongestionSolver]
    ],
    First
];

(* --- Private Section --- *)
Begin["`Private`"];

(* Defaults *)
$MFGraphsVerbose = False;
$MFGraphsParallelThreshold = 6;
$MFGraphsParallelReady = False;
$MFGraphsParallelLaunchThreshold = 50;

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
MFGPrint[args___] := If[$MFGraphsVerbose, Print[args]];
MFGPrintTemporary[args___] := If[$MFGraphsVerbose, PrintTemporary[args], Null];
(* :!CodeAnalysis::EndBlock:: *)

EnsureParallelKernels[] := If[$KernelCount === 0, LaunchKernels[]];

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

(* --- Shared infrastructure for FictitiousPlay and Monotone methods --- *)

Hess[j_?NumberVectorQ] := DiagonalMatrix[1/# & /@ j];

InverseHessian[j_?NumberVectorQ] := DiagonalMatrix[j];

LookupAssociationValue[assoc_Association, key_, default_:0.] :=
  If[KeyExistsQ[assoc, key], assoc[key], default];

BuildMonotoneStateData[d2e_Association] :=
  Module[{B, KM, jj, halfPairs, signedFlows, rules, substituted, S},
    {B, KM, jj} = GetKirchhoffLinearSystem[d2e];
    halfPairs = List @@@ Lookup[d2e, "edgeList", {}];
    If[Length[jj] === 0 || Length[halfPairs] === 0,
      S = SparseArray[{}, {Length[halfPairs], Length[jj]}],
      signedFlows = Lookup[d2e, "SignedFlows", <||>];
      rules = Join[
        Lookup[d2e, "RuleBalanceGatheringFlows", <||>],
        Lookup[d2e, "RuleExitFlowsIn", <||>],
        Lookup[d2e, "RuleEntryOut", <||>]
      ];
      substituted = (signedFlows[#] & /@ halfPairs) /. rules;
      S = Last @ CoefficientArrays[substituted, jj]
    ];
    <|
      "B" -> B,
      "KM" -> KM,
      "FullVariables" -> jj,
      "HalfPairs" -> halfPairs,
      "SignedEdgeMatrix" -> S
    |>
  ];

BuildMonotoneValueSystem[d2e_Association] :=
  Module[{statePairs, pairIndex, equalityPairs, boundaryValues, transitions, switching,
      stateCount, transitionsAt},
    statePairs = DeleteDuplicates @ Cases[Lookup[d2e, "us", {}], u[a_, b_] :> {a, b}];
    stateCount = Length[statePairs];
    If[stateCount === 0,
      Return[<|"StateValueAssociation" -> (<||> &)|>, Module]
    ];
    pairIndex = AssociationThread[statePairs, Range[stateCount]];
    equalityPairs = Cases[Lookup[d2e, "Equ", {}], (u @@@ p1_) == (u @@@ p2_) :> {p1, p2}];
    boundaryValues = Cases[Lookup[d2e, "Equ", {}], (u @@@ p_) == val_?NumericQ :> {p, val}];
    transitions = Lookup[d2e, "Transitions", <||>];
    switching = Lookup[d2e, "SwitchingCosts", <||>];
    transitionsAt = GroupBy[Lookup[d2e, "auxTriples", {}], #[[{1, 2}]] &];

    <|
      "StateValueAssociation" -> Function[{edgeCosts},
        Module[{vars, eqs, sol, usLocal},
          vars = u @@@ statePairs;
          eqs = Flatten @ {
            Map[
              Function[{pair},
                u @@ pair == (
                  LookupAssociationValue[transitions, pair, 0] +
                  Plus @@ Map[
                    Function[{triple},
                      u @@ triple[[{2, 3}]] + LookupAssociationValue[switching, triple, 0]
                    ],
                    Lookup[transitionsAt, Key[pair], {}]
                  ]
                )
              ],
              statePairs
            ],
            Map[(u @@@ #[[1]] == u @@@ #[[2]]) &, equalityPairs],
            Map[(u @@@ #[[1]] == #[[2]]) &, boundaryValues]
          };
          usLocal = Lookup[d2e, "us", {}];
          sol = Quiet @ Check[Solve[eqs, usLocal], {}];
          If[MatchQ[sol, {{_Rule ..}}],
            Association[First[sol]],
            Missing["NotAvailable"]
          ]
        ]
      ]
    |>
  ];

BuildReducedKirchhoffCoordinates[d2e_Association, basePoint_:Automatic] :=
  Module[{stateData, KM, jj, S, nullBasis, invisibleBasis, visibleBasis,
      removeInvisible, basisMatrix, dim, baseVec, edgeOffset, edgeBasis,
      invisibleDim, signedEdgeRows, tol = 10^-10},
    stateData = BuildMonotoneStateData[d2e];
    KM = stateData["KM"];
    jj = stateData["FullVariables"];
    S = stateData["SignedEdgeMatrix"];
    nullBasis = N @ Orthogonalize[NullSpace[Normal[KM]]];
    invisibleBasis = N @ Orthogonalize[NullSpace[Join[Normal[KM], Normal[S]]]];
    removeInvisible[v_] :=
      Fold[#1 - (#1.#2) #2 &, v, invisibleBasis];
    visibleBasis = If[nullBasis === {},
      {},
      If[invisibleBasis === {},
        nullBasis,
        Orthogonalize @ Select[removeInvisible /@ nullBasis, Norm[#] > tol &]
      ]
    ];
    visibleBasis = Select[visibleBasis, Norm[#] > tol &];
    basisMatrix = If[visibleBasis === {},
      ConstantArray[0., {Length[jj], 0}],
      N @ Transpose[visibleBasis]
    ];
    dim = Length[visibleBasis];
    baseVec = Which[
      basePoint === Automatic && Length[jj] === 0, {},
      basePoint === Automatic, Quiet @ Check[
        N @ LeastSquares[N @ Normal[KM], N @ stateData["B"]],
        Missing["NotAvailable"]
      ],
      True, N @ basePoint
    ];
    signedEdgeRows = If[ArrayQ[S], First @ Dimensions[S], 0];
    edgeOffset = If[signedEdgeRows === 0,
      {},
      If[ListQ[baseVec] && VectorQ[baseVec, NumericQ],
        N @ (S . baseVec),
        Missing["NotAvailable"]
      ]
    ];
    edgeBasis = If[signedEdgeRows === 0,
      ConstantArray[0., {0, dim}],
      If[MatrixQ[basisMatrix], N @ (S . basisMatrix), Missing["NotAvailable"]]
    ];
    invisibleDim = Length[invisibleBasis];
    <|
      "BasePoint" -> baseVec,
      "BasisMatrix" -> basisMatrix,
      "StateDimension" -> dim,
      "FullDimension" -> Length[jj],
      "CostInvisibleDimension" -> invisibleDim,
      "SignedEdgeOffset" -> edgeOffset,
      "SignedEdgeBasis" -> edgeBasis,
      "FullVariables" -> jj
    |>
  ];

BuildMonotonePairCostAssociation[halfPairs_List, edgeList_List, q_] :=
  Association @ Flatten @ MapThread[
    Function[{pair, edge, flow},
      With[{cost = If[PossibleZeroQ[flow], 0., N @ (Abs[flow]^(1) + 0)]}, (* Critical alpha=1 *)
        {pair -> cost, Reverse[pair] -> cost}
      ]
    ],
    {halfPairs, edgeList, q}
  ];

BuildCriticalQuadraticEdgeModel[d2e_Association, tol_:10^-8] :=
  Module[{edges},
    edges = Lookup[d2e, "edgeList", {}];
    (* Critical congestion always assumes alpha=1, so slopes are all 1 *)
    AssociationThread[edges, ConstantArray[1., Length[edges]]]
  ];

BuildCriticalQuadraticObjective[d2e_Association] :=
  Module[{stateData, B, KM, jj, S, q, edgeSlopes, switching, exitRules,
      outExitPairs, exitCostByPair, linearTerm, quadraticTerm, edgeSlopeVector},
    stateData = BuildMonotoneStateData[d2e];
    B = stateData["B"];
    KM = stateData["KM"];
    jj = stateData["FullVariables"];
    If[Length[jj] === 0,
      Return[$Failed, Module]
    ];
    edgeSlopes = BuildCriticalQuadraticEdgeModel[d2e];
    If[edgeSlopes === $Failed,
      Return[$Failed, Module]
    ];
    S = stateData["SignedEdgeMatrix"];
    q = If[Length[stateData["HalfPairs"]] === 0, {}, S . jj];
    edgeSlopeVector = Lookup[edgeSlopes, Lookup[d2e, "edgeList", {}], 1.];
    quadraticTerm =
      If[q === {} || edgeSlopeVector === {},
        0,
        1/2 q . DiagonalMatrix[edgeSlopeVector] . q
      ];
    switching = Lookup[d2e, "SwitchingCosts", <||>];
    exitRules = Lookup[d2e, "RuleExitValues", <||>];
    outExitPairs = List @@@ Lookup[d2e, "exitEdges", {}];
    exitCostByPair =
      Association @ Table[
        pair -> LookupAssociationValue[exitRules, u @@ pair, 0],
        {pair, Join[outExitPairs, Reverse /@ outExitPairs]}
      ];
    linearTerm =
      Total[
        (
          LookupAssociationValue[switching, #, 0] +
          LookupAssociationValue[exitCostByPair, #[[{2, 3}]], 0]
        ) (j @@ #) & /@ Lookup[d2e, "auxTriples", {}]
      ] +
      10^-7 Total[jj];
    <|
      "B" -> B,
      "KM" -> KM,
      "Variables" -> jj,
      "Objective" -> quadraticTerm + linearTerm,
      "StateData" -> stateData
    |>
  ];

UseQuadraticCriticalBackendQ[d2e_Association] :=
  Module[{edges, graph, degrees, maxDegree, nonzeroSwitching},
    edges = Lookup[d2e, "edgeList", {}];
    graph = Graph[edges];
    degrees = VertexDegree[graph];
    maxDegree = If[degrees === {}, 0, Max[degrees]];
    nonzeroSwitching =
      AnyTrue[
        Values @ Lookup[d2e, "SwitchingCosts", <||>],
        NumericQ[#] && !PossibleZeroQ[#] &
      ];
    TrueQ[AcyclicGraphQ[graph]] && (maxDegree <= 2 || !nonzeroSwitching)
  ];

MonotoneVariableFieldValue[var_, values_Association, switching_Association] :=
  Replace[
    var,
    {
      j[r_, i_, w_] :> N @ (
        LookupAssociationValue[switching, {r, i, w}, 0] +
        LookupAssociationValue[values, u[i, w], 0.] +
        LookupAssociationValue[values["PairCosts"], {i, w}, 0.]
      ),
      j[a_, b_] :> N @ (
        LookupAssociationValue[values, u[a, b], 0.] +
        LookupAssociationValue[values["PairCosts"], {a, b}, 0.]
      )
    }
  ];

EncodeFlowAssociation[numericState_Association, assoc_Association] :=
    Module[{vars, vals},
        vars = Lookup[numericState, "FlowVariables", {}];
        If[vars === {},
            Return[Developer`ToPackedArray[{}], Module]
        ];
        vals = vars /. assoc;
        (* Replace any remaining symbolic MFG variables with 0.0 *)
        vals = vals /. {u[___] :> 0., j[___] :> 0., z[___] :> 0.};
        vals = Quiet @ Check[N @ vals, vals];
        If[VectorQ[vals, NumericQ],
            Developer`ToPackedArray @ N @ vals,
            Missing["NotAvailable"]
        ]
    ];

EncodeFlowAssociation[_, _] := Missing["NotAvailable"];

DecodeFlowVector[numericState_Association, flowVec_] :=
    Module[{vars, vals},
        vars = Lookup[numericState, "FlowVariables", {}];
        vals = Quiet @ Check[Developer`ToPackedArray @ N @ flowVec, $Failed];
        If[vals === $Failed || !VectorQ[vals, NumericQ] || Length[vals] =!= Length[vars],
            <||>,
            AssociationThread[vars, vals]
        ]
    ];

DecodeFlowVector[_, _] := <||>;

ComputeSignedEdgeFlowsFast[numericState_Association, flowVec_] :=
    Module[{indices, signedEdgeMatrix, kirchhoffVec},
        indices = Lookup[numericState, "KirchhoffIndicesInFlowVector", {}];
        signedEdgeMatrix = Lookup[numericState, "SignedEdgeMatrix", SparseArray[{}, {0, 0}]];
        If[Dimensions[signedEdgeMatrix][[1]] === 0,
            Return[Developer`ToPackedArray[{}], Module]
        ];
        If[indices === {} || !VectorQ[indices, IntegerQ],
            Return[Missing["NotAvailable"], Module]
        ];
        If[Length[flowVec] < Max[indices],
            Return[Missing["NotAvailable"], Module]
        ];
        kirchhoffVec = flowVec[[indices]];
        If[!VectorQ[kirchhoffVec, NumericQ],
            Return[Missing["NotAvailable"], Module]
        ];
        Developer`ToPackedArray @ N @ (signedEdgeMatrix . kirchhoffVec)
    ];

ComputeKirchhoffResidualFast[numericState_Association, flowVec_] :=
    Module[{indices, KM, B, kirchhoffVec, residual},
        indices = Lookup[numericState, "KirchhoffIndicesInFlowVector", {}];
        KM = Lookup[numericState, "KirchhoffKM", SparseArray[{}, {0, 0}]];
        B = Lookup[numericState, "KirchhoffB", Developer`ToPackedArray[{}]];
        If[indices === {} || !VectorQ[indices, IntegerQ],
            Return[Missing["NotAvailable"], Module]
        ];
        If[Length[flowVec] < Max[indices],
            Return[Missing["NotAvailable"], Module]
        ];
        kirchhoffVec = flowVec[[indices]];
        If[!VectorQ[kirchhoffVec, NumericQ],
            Return[Missing["NotAvailable"], Module]
        ];
        residual = Quiet @ Check[N @ Norm[KM . kirchhoffVec - B, Infinity], Missing["NotAvailable"]];
        If[NumericQ[residual], residual, Missing["NotAvailable"]]
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
                (j @@@ entryPairs) /. flowAssoc
            ];
        exitFlows =
            If[exitPairs === {},
                {},
                (j @@@ exitPairs) /. flowAssoc
            ];
        (* Default symbolic variables to 0.0 for mass calculation *)
        entryFlows = entryFlows /. {j[___] :> 0., u[___] :> 0.};
        exitFlows = exitFlows /. {j[___] :> 0., u[___] :> 0.};
        
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

RecoverCriticalFlowAssociation[Eqs_Association, flowBase_Association, tol_?NumericQ] :=
    Module[{entryInRules, flowRules, replacementRules, allFlowVars, resolve, resolved},
        entryInRules = Association @ Flatten[ToRules /@ Lookup[Eqs, "EqEntryIn", {}]];
        flowRules = Join[
            entryInRules,
            Lookup[Eqs, "RuleEntryOut", <||>],
            Lookup[Eqs, "RuleExitFlowsIn", <||>],
            Lookup[Eqs, "RuleBalanceGatheringFlows", <||>]
        ];
        replacementRules = Join[Normal[flowBase], Normal[flowRules]];
        resolve[expr_] := FixedPoint[
            Function[val, Quiet @ Check[val /. replacementRules, val]],
            expr,
            10
        ];
        allFlowVars = Join[Lookup[Eqs, "js", {}], Lookup[Eqs, "jts", {}]];
        If[!ListQ[allFlowVars] || allFlowVars === {},
            Return[<||>, Module]
        ];
        resolved = Association @ Cases[
            Map[
                Function[{var},
                    Module[{val},
                        val = Quiet @ Check[N @ resolve[var], Missing["NotAvailable"]];
                        If[NumericQ[val],
                            var -> If[Abs[val] <= tol, 0., N[val]],
                            Nothing
                        ]
                    ]
                ],
                allFlowVars
            ],
            _Rule
        ];
        If[
            !AssociationQ[resolved] ||
            !And @@ (KeyExistsQ[resolved, #] & /@ allFlowVars),
            Failure["RecoverFlow", <|"Reason" -> "FlowReconstructionFailed"|>],
            resolved
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
            VectorQ[bEqRaw, NumericQ], bEqRaw,
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
            "KirchhoffResidual" -> residual,
            "NonnegativityViolation" -> nnViolation,
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
            "UtilityAssociation" -> utilityAssoc,
            "NodePotentials" -> nodePotentials,
            "ReducedCostAssociation" -> reducedCostAssoc,
            "EdgeCostAssociation" -> edgeCostAssoc,
            "BellmanResidual" -> bellmanResidual,
            "CostNonnegativeQ" -> True,
            "Residuals" -> <|
                "BellmanResidual" -> bellmanResidual
            |>
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

        topoData = Lookup[staticData, "TopologicalOrder", Missing["NotAvailable"]];
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
            VectorQ[bEqRaw, NumericQ], bEqRaw,
            True, {}
        ];
        massVec = Developer`ToPackedArray @ N @ Lookup[newFlowAssoc, massVars, 0.];
        residual = Quiet @ Check[N @ Norm[aEq . massVec - bEq, Infinity], Infinity];
        nnViolation = Quiet @ Check[N @ Max[0., -Min[newFlowVec]], Infinity];

        newNodeMassAssoc = AssociationThread[statePairs, Lookup[newFlowAssoc, jVars, 0.]];

        edgeListPairs = List @@@ Lookup[eqs, "edgeList", {}];
        signedVals = Lookup[newFlowAssoc, j @@@ edgeListPairs, 0.];
        comparableFlowVec = signedVals;

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
                "SpatialFlowAssociation" -> <||>,
                "TransitionFlowAssociation" -> <||>,
                "NodeMassAssociation" -> newNodeMassAssoc,
                "ComparableFlowVector" -> comparableFlowVec,
                "FeasibleQ" -> True,
                "KirchhoffResidual" -> residual,
                "NonnegativityViolation" -> nnViolation,
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
            "SupportSignature" -> currentSignature,
            "StableIterations" -> currentStableCount,
            "StableQ" -> stableQ,
            "OracleReadyQ" -> stableQ,
            "ActiveVariables" -> activeVars,
            "InactiveVariables" -> inactiveVars,
            "AmbiguousVariables" -> ambiguousVars,
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
        <|"Convergence" -> Missing["NotApplicable"]|>,
        extraAssoc
    ];

SolveMFGCompiledInputQ[input_Association] :=
    KeyExistsQ[input, "EqGeneral"] &&
    KeyExistsQ[input, "js"] &&
    KeyExistsQ[input, "jts"];

SolveMFGCompileInput[input_Association, opts_List : {}] :=
    If[SolveMFGCompiledInputQ[input],
        input,
        MFGraphs`DataToEquations[input]
    ];

CheckFlowFeasibility[Null] := "Infeasible";
CheckFlowFeasibility[assoc_Association] :=
    Module[{flowKeys, flowVals},
        flowKeys = Select[Keys[assoc], MatchQ[#, _j] &];
        flowVals = Lookup[assoc, flowKeys];
        If[flowVals === {} || Min[Select[flowVals, NumericQ]] < 0,
            "Infeasible", "Feasible"]
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

SolveMFGRunStage[
    method_String,
    input_Association,
    opts_List,
    dispatchMode_String : "Direct"
] :=
    Module[{d2e},
        Switch[method,
            "CriticalCongestion",
                d2e = If[SolveMFGCompiledInputQ[input], input, SolveMFGCompileInput[input, opts]];
                CriticalCongestionSolver[
                    d2e,
                    Sequence @@ FilterRules[opts, Options[CriticalCongestionSolver]]
                ]
            ,
            _,
                SolveMFGUnknownMethodResult[method]
        ]
    ];

SolveMFGAutomaticDispatch[input_Association, opts_List] :=
    Module[{d2e, result = Missing["NotAvailable"], decision, trace = {}, methodUsed = "Automatic",
            attempt = "CriticalCongestion"},
        d2e = If[SolveMFGCompiledInputQ[input], input, Quiet @ SolveMFGCompileInput[input, opts]];
        If[!AssociationQ[d2e],
            Return[SolveMFGAbortResult["DataToEquationsFailed", methodUsed, trace], Module]
        ];
        result = SolveMFGRunStage[attempt, d2e, opts, "Automatic"];
        decision = SolveMFGClassifyOutcome[result];
        methodUsed = attempt;
        trace = Append[trace, SolveMFGBuildTraceEntry[methodUsed, result, decision]];
        Switch[decision,
            "Done" | "Fallback",
                If[AssociationQ[result],
                    SolveMFGFinalizeWithTrace[result, methodUsed, trace],
                    SolveMFGAbortResult["InvalidSolverReturnShape", methodUsed, trace]
                ],
            _,
                SolveMFGAbortResult["InvalidSolverReturnShape", methodUsed, trace]
        ]
    ];

SolveMFG[input_Association, opts:OptionsPattern[]] :=
    Module[{method, normalizedMethod},
        method = OptionValue[Method];
        normalizedMethod = ToString[method];

        Switch[normalizedMethod,
            "CriticalCongestion" | "CriticalCongestionSolver",
                SolveMFGRunStage["CriticalCongestion", input, {opts}, "Direct"]
            ,
            "Automatic",
                SolveMFGAutomaticDispatch[input, {opts}]
            ,
            _,
                SolveMFGUnknownMethodResult[method]
        ]
    ];

SolveMFG[_, opts:OptionsPattern[]] :=
    SolveMFGUnknownMethodResult[OptionValue[Method]];

End[];

EndPackage[];
