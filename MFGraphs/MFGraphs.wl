(* Wolfram Language Package *)
(*
   MFGraphs: A Wolfram Language package for Mean Field Games on Networks.

   This top-level file is a thin bootstrap.  It owns the BeginPackage/
   EndPackage block, the verbosity/parallel flags, a handful of shared
   public helpers (MFGPrint, MFGParallelMap, ...), and a small private
   library of result-envelope builders used by every solver.  All solver
   logic lives in the submodules loaded at the bottom of this file.

   Main Components:
   - Examples/*            - built-in test cases.
   - DNFReduce             - symbolic logical reduction engine.
   - DataToEquations       - symbolic graph-to-equation converter.
   - Solvers               - critical-congestion solver suite.
   - NonLinearSolver       - iterative solver for general congestion.
   - Monotone              - Hessian Riemannian Flow numerical solver.
   - TimeDependentSolver   - backward-forward sweep solver.
   - Graphics              - public visualization helpers.
   - Scenario              - typed scenario kernel.
   - FictitiousPlayBackend - Phase 5 numeric backend (internal).
   - SolveMFGDispatch      - unified SolveMFG entrypoint and cascade.
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
or Length[list] >= $MFGraphsParallelLaunchThreshold (when kernels must be launched).
This avoids paying ~3s kernel launch overhead for small workloads.";

$MFGraphsParallelLaunchThreshold::usage =
"$MFGraphsParallelLaunchThreshold is the minimum list length required to justify
launching parallel subkernels. Only applies when $KernelCount === 0. Default is 50.";

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

(* Load submodules in dependency order *)
Get["MFGraphs`DNFReduce`"];
Get["MFGraphs`Scenario`"];
Get["MFGraphs`Examples`ExampleScenarios`"];
Get["MFGraphs`Examples`ExamplesData`"];
Get["MFGraphs`Unknowns`"];
Get["MFGraphs`System`"];
Get["MFGraphs`DataToEquations`"];
Get["MFGraphs`Solvers`"];
Get["MFGraphs`Graphics`"];

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
            {B, KM, jj} = GetKirchhoffLinearSystem[Eqs];
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
        <|"Convergence" -> Missing["NotApplicable"]|>,
        extraAssoc
    ];

SolveMFGCompiledInputQ[input_Association] :=
    KeyExistsQ[input, "EqGeneral"] &&
    KeyExistsQ[input, "js"] &&
    KeyExistsQ[input, "jts"];

(* Load submodules in dependency order *)
Scan[
    Get,
    {
        "MFGraphs`Examples`ExamplesData`",
        "MFGraphs`DNFReduce`",
        "MFGraphs`DataToEquations`",
        "MFGraphs`Solvers`",
        "MFGraphs`Graphics`",
        "MFGraphs`FictitiousPlayBackend`",
        "MFGraphs`SolveMFGDispatch`"
    }
];

(* Soft-optional legacy modules: load only when present on $Path. *)
Scan[
    Function[ctx,
        If[StringQ[Quiet @ Check[FindFile[ctx], $Failed]],
            Get[ctx]
        ]
    ],
    {
        "MFGraphs`NonLinearSolver`",
        "MFGraphs`Monotone`",
        "MFGraphs`TimeDependentSolver`",
        "MFGraphs`Examples`TimeDependentExamples`"
    }
];

End[];

EndPackage[];
