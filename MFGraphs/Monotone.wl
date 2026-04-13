(* Wolfram Language package *)
(*
   Monotone: Numerical solver for Multi-Population Wardrop Equilibrium.
   
   Theoretical Basis:
   - "Hessian Riemannian Flow For Multi-Population Wardrop Equilibrium" (Bakaryan et al., 2025)
   
   This module implements an efficient numerical solver using the Hessian Riemannian Flow (HRF) 
   method. It reformulates the equilibrium problem as a distributed optimization problem 
   projected onto the manifold defined by flow conservation constraints.
*)

(* --- Public API declarations --- *)

Hess::usage =
"Hess[j] returns a numeric diagonal matrix with the reciprocals of the elements in the (numeric) vector j.";

InverseHessian::usage =
"InverseHessian[j] returns a diagonal matrix with the (numeric) vector j. This is the inverse of the matrix Hess[j].";

NumberMatrixQ::usage =
"NumberMatrixQ[A] returns True if the elements of the matrix A are numeric.";

HessianSandwich::usage =
"HessianSandwich[j, A, At] returns the product A . InverseHessian[j] . At.";

GradientProjection::usage =
"GradientProjection[x, A, dim, At] returns the projected gradient operator matrix.
This is InverseHessian[x] . (I - At . PseudoInverse[A . InverseHessian[x] . At] . A . InverseHessian[x]).";

CachedGradientProjection::usage =
"CachedGradientProjection[x, KM, dim, At, cache] is a version of GradientProjection
that caches the projected linear solve and reuses it when x has not changed significantly.
cache must be a symbol holding a 1-element list {Null} or
{<|\"x\" -> ..., \"Method\" -> ..., \"Solver\" -> ...|>}.
The function has the HoldAll attribute so that cache is passed by reference.";

MonotoneSolverFromData::usage =
"MonotoneSolverFromData[Data] solves the MFG problem from raw Data using the monotone operator method.
Options: \"ResidualTolerance\" (default 10^-6), \"MaxTime\" (default 100), \
\"MaxSteps\" (default 5000), \"UseCachedProjection\" (default True), \
\"PotentialFunction\", \"CongestionExponentFunction\", and \"InteractionFunction\" \
(default Automatic, meaning use the current global MFGraphs` definitions of V, alpha, and g). \
It always returns a standardized solver-result association containing solver metadata, \
comparison fields, convergence data, and the solver-specific payload key \"AssoMonotone\".";

MonotoneSolver::usage =
"MonotoneSolver[d2e] solves the MFG problem using the monotone operator method.
Options: \"ResidualTolerance\" (default 10^-6), \"MaxTime\" (default 100), \
\"MaxSteps\" (default 5000), \"UseCachedProjection\" (default True), \
\"PotentialFunction\", \"CongestionExponentFunction\", and \"InteractionFunction\" \
(default Automatic, meaning use the current global MFGraphs` definitions of V, alpha, and g). \
It always returns a standardized solver-result association containing solver metadata, \
comparison fields, convergence data, and the solver-specific payload key \"AssoMonotone\".";

MonotoneSolverODE::usage =
"MonotoneSolverODE[reducedState, KM, B, cc] solves the reduced Hessian Riemannian
flow and returns an association with the reconstructed solution and convergence data.
Options: \"ResidualTolerance\" (default 10^-6), \"MaxTime\" (default 100),
\"MaxSteps\" (default 5000), and \"UseCachedProjection\" (default True).";

MonotoneSolver::degenerate =
"Degenerate case: no flow variables were found, so the monotone solver was skipped.";

MonotoneSolver::seedfail =
"Failed to find an interior feasible seed for the monotone ODE solve.";

MonotoneSolver::odefail =
"The monotone ODE solve failed before a usable partial trajectory was produced.";

Options[MonotoneSolverFromData] = {
    "ResidualTolerance" -> 10^-6,
    "MaxTime" -> 100,
    "MaxSteps" -> 5000,
    "UseCachedProjection" -> True,
    "PotentialFunction" -> Automatic,
    "CongestionExponentFunction" -> Automatic,
    "InteractionFunction" -> Automatic
};
Options[MonotoneSolver] = Options[MonotoneSolverFromData];
Options[MonotoneSolverODE] = {
    "ResidualTolerance" -> 10^-6,
    "MaxTime" -> 100,
    "MaxSteps" -> 5000,
    "UseCachedProjection" -> True
};

Begin["`Private`"];

(* HoldAll so that the cache argument (a mutable list) is passed by reference *)
SetAttributes[CachedGradientProjection, HoldAll];

(* --- Matrix utilities --- *)

Hess[j_?NumberVectorQ] := DiagonalMatrix[1/# & /@ j];

InverseHessian[j_?NumberVectorQ] := DiagonalMatrix[j];

NumberMatrixQ[A_] := NumberVectorQ@Flatten@A;

HessianSandwich[x_?NumberVectorQ, A_, At_] := A . InverseHessian[x] . At;
HessianSandwich[x_?NumberVectorQ, A_] := A . InverseHessian[x] . Transpose[A];

BuildProjectionCacheEntry[x_?NumberVectorQ, A_, At_] :=
  Module[{invH, sandwich, solver},
    If[Length[x] === 0,
      Return[<|"x" -> x, "Method" -> "Degenerate"|>, Module]
    ];
    invH = InverseHessian[x];
    sandwich = HessianSandwich[x, A, At];
    solver = Quiet @ Check[LinearSolve[sandwich], $Failed];
    If[solver === $Failed,
      <|"x" -> x, "Method" -> "PseudoInverse", "PseudoInverse" -> PseudoInverse[sandwich]|>,
      <|"x" -> x, "Method" -> "LinearSolve", "Solver" -> solver|>
    ]
  ];

ProjectionOperatorFromEntry[x_?NumberVectorQ, A_, dim_Integer, At_, entry_Association] :=
  If[dim === 0,
    {},
    Module[{invH, correction},
    invH = InverseHessian[x];
    correction = If[Lookup[entry, "Method", "PseudoInverse"] === "LinearSolve",
      At . entry["Solver"][A . invH],
      At . entry["PseudoInverse"] . A . invH
    ];
    invH . (IdentityMatrix[dim] - correction)
    ]
  ];

(* --- GradientProjection: projected gradient operator --- *)

GradientProjection[x_?NumberVectorQ, A_, dim_, At_] :=
  ProjectionOperatorFromEntry[x, A, dim, At, BuildProjectionCacheEntry[x, A, At]];

GradientProjection[x_?NumberVectorQ, A_] :=
	Module[{dim, At},
		dim = Length[x];
		At = Transpose[A];
		ProjectionOperatorFromEntry[x, A, dim, At, BuildProjectionCacheEntry[x, A, At]]
	];

(* --- CachedGradientProjection: caches PseudoInverse when x changes slowly --- *)
(* HoldAll keeps cache as the caller's symbol so part-assignment works. *)
(* Other arguments are force-evaluated via Module initializers. *)

CachedGradientProjection[x_, KM_, dim_, At_, cache_, tol_:10^-6] :=
  Module[{xVal = x, KMVal = KM, dimVal = dim, AtVal = At, tolVal = tol, entry},
    If[cache[[1]] =!= Null && Norm[xVal - cache[[1]]["x"], Infinity] < tolVal,
      entry = cache[[1]],
      entry = BuildProjectionCacheEntry[xVal, KMVal, AtVal];
      cache[[1]] = entry
    ];
    ProjectionOperatorFromEntry[xVal, KMVal, dimVal, AtVal, entry]
  ];

(* --- BuildMonotoneStateData: shared linear state for MonotoneSolver --- *)

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

ConjunctionTerms[expr_] :=
  Which[
    TrueQ[expr], {},
    Head[expr] === And, List @@ expr,
    True, {expr}
  ];

BuildMonotoneValueSystem[d2e_Association] :=
  Module[{statePairs, pairIndex, equalityPairs, boundaryValues, transitions, switching,
      stateCount},
    statePairs = DeleteDuplicates @ Cases[Lookup[d2e, "us", {}], u[a_, b_] :> {a, b}];
    stateCount = Length[statePairs];
    If[stateCount === 0,
      Return[<|"StateValueAssociation" -> (<||> &)|>, Module]
    ];
    pairIndex = AssociationThread[statePairs, Range[stateCount]];
    equalityPairs = DeleteDuplicates @ Select[
      Cases[
        Join[
          KeyValueMap[Equal, Lookup[d2e, "RuleEntryValues", <||>]],
          ConjunctionTerms @ Lookup[d2e, "EqValueAuxiliaryEdges", True]
        ],
        Equal[u[a_, b_], u[c_, d_]] :> {{a, b}, {c, d}}
      ],
      KeyExistsQ[pairIndex, #[[1]]] && KeyExistsQ[pairIndex, #[[2]]] &
    ];
    boundaryValues = Association @ KeyValueMap[
      (#1 -> N[#2]) &,
      Lookup[d2e, "RuleExitValues", <||>]
    ];
    transitions = Select[
      Lookup[d2e, "auxTriples", {}],
      KeyExistsQ[pairIndex, #[[{1, 2}]]] && KeyExistsQ[pairIndex, #[[{2, 3}]]] &
    ];
    switching = Lookup[d2e, "SwitchingCosts", <||>];
    <|
      "StateValueAssociation" ->
        Function[{pairCosts},
          Module[{dist, changed, sourceIdx, targetIdx, newValue, iter = 0},
            dist = ConstantArray[Infinity, stateCount];
            KeyValueMap[
              Function[{key, value},
                If[KeyExistsQ[pairIndex, List @@ key],
                  dist[[pairIndex[List @@ key]]] = Min[dist[[pairIndex[List @@ key]]], value]
                ]
              ],
              boundaryValues
            ];
            While[iter < stateCount,
              iter++;
              changed = False;
              Do[
                sourceIdx = pairIndex[equality[[1]]];
                targetIdx = pairIndex[equality[[2]]];
                If[dist[[targetIdx]] < dist[[sourceIdx]],
                  dist[[sourceIdx]] = dist[[targetIdx]];
                  changed = True
                ];
                If[dist[[sourceIdx]] < dist[[targetIdx]],
                  dist[[targetIdx]] = dist[[sourceIdx]];
                  changed = True
                ],
                {equality, equalityPairs}
              ];
              Do[
                sourceIdx = pairIndex[triple[[{1, 2}]]];
                targetIdx = pairIndex[triple[[{2, 3}]]];
                newValue =
                  LookupAssociationValue[switching, triple, 0.] +
                  LookupAssociationValue[pairCosts, triple[[{2, 3}]], 0.] +
                  dist[[targetIdx]];
                If[newValue < dist[[sourceIdx]],
                  dist[[sourceIdx]] = newValue;
                  changed = True
                ],
                {triple, transitions}
              ];
              If[!changed, Break[]]
            ];
            AssociationThread[u @@@ statePairs, N @ dist]
          ]
        ]
    |>
  ];

LookupAssociationValue[assoc_Association, key_, default_:0.] :=
  If[KeyExistsQ[assoc, key], assoc[key], default];

BuildCriticalQuadraticEdgeModel[d2e_Association, tol_:10^-8] :=
  Module[{edges, alphaValues, samples, slopes},
    edges = Lookup[d2e, "edgeList", {}];
    alphaValues = Quiet @ Check[N[alpha /@ edges], $Failed];
    If[
      alphaValues =!= $Failed &&
      VectorQ[alphaValues, NumericQ] &&
      Max[Abs[alphaValues - 1.]] <= tol,
      Return[AssociationThread[edges, ConstantArray[1., Length[edges]]], Module]
    ];
    samples =
      Quiet @ Check[
        Table[
          {N @ Cost[1., edge], N @ Cost[2., edge], N @ Cost[3., edge]},
          {edge, edges}
        ],
        $Failed
      ];
    If[samples === $Failed || !MatrixQ[samples, NumericQ],
      Return[$Failed, Module]
    ];
    slopes = Map[
      Function[{vals},
        Module[{c1, c2, c3, slope, intercept},
          {c1, c2, c3} = vals;
          slope = c2 - c1;
          intercept = 2 c1 - c2;
          If[
            slope < -tol ||
            Abs[intercept] > tol ||
            Abs[(3 slope + intercept) - c3] > 10 tol,
            Return[$Failed, Module]
          ];
          N @ slope
        ]
      ],
      samples
    ];
    If[MemberQ[slopes, $Failed],
      $Failed,
      AssociationThread[edges, slopes]
    ]
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

BuildMonotoneComparisonData[d2e_Association, stateData_Association, solution_Association] :=
  Module[{missing, edgeList, flowAssoc, jj, jjVals, signedVals, residual,
      boundaryMassData, utilityReductionResidualData},
    missing = Missing["NotAvailable"];
    edgeList = Lookup[d2e, "edgeList", {}];
    flowAssoc = Association @ KeySelect[solution, MatchQ[#, _j] &];
    boundaryMassData = BuildBoundaryMassData[d2e, flowAssoc];
    utilityReductionResidualData = BuildUtilityReductionResidualData[d2e, solution];
    jj = stateData["FullVariables"];
    jjVals = Lookup[flowAssoc, jj, missing];
    signedVals =
      If[ListQ[jjVals] && VectorQ[jjVals, NumericQ],
        If[Length[edgeList] === 0, {}, N @ (stateData["SignedEdgeMatrix"] . jjVals)],
        missing
      ];
    residual =
      If[ListQ[jjVals] && VectorQ[jjVals, NumericQ],
        N @ Norm[stateData["KM"] . jjVals - stateData["B"], Infinity],
        missing
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

BuildReducedKirchhoffMetadata[stateData_Association] :=
  Module[{KM, S, fullDim, kmRank, stackedRank},
    KM = Normal @ Lookup[stateData, "KM", {}];
    S = Normal @ Lookup[stateData, "SignedEdgeMatrix", {}];
    fullDim = Length @ Lookup[stateData, "FullVariables", {}];
    If[fullDim === 0,
      Return[
        <|
          "ReducedStateDimension" -> 0,
          "FullStateDimension" -> 0,
          "CostInvisibleDimension" -> 0
        |>,
        Module
      ]
    ];
    kmRank = MatrixRank[KM];
    stackedRank = MatrixRank[Join[KM, S]];
    <|
      "ReducedStateDimension" -> Max[0, stackedRank - kmRank],
      "FullStateDimension" -> fullDim,
      "CostInvisibleDimension" -> Max[0, fullDim - stackedRank]
    |>
  ];

CleanCriticalQuadraticSolution[model_Association, rawSolution_List] :=
  Module[{vars, S, qTarget, cleaned},
    vars = model["Variables"];
    S = model["StateData"]["SignedEdgeMatrix"];
    qTarget = If[
      Length[Lookup[model["StateData"], "HalfPairs", {}]] === 0,
      {},
      N @ (S . (vars /. rawSolution))
    ];
    cleaned =
      Quiet @ Check[
        LinearOptimization[
          Total[vars],
          Join[
            Thread[Normal[model["KM"]] . vars == Normal[model["B"]]],
            If[qTarget === {}, {}, Thread[S . vars == qTarget]],
            Thread[vars >= 0]
          ],
          vars
        ],
        $Failed
      ];
    If[MatchQ[cleaned, {_Rule ..}], cleaned, rawSolution]
  ];

TryCriticalQuadraticMonotoneSolve[d2e_Association, residualTolerance_:10^-6] :=
  Module[{model, vars, rawSolution, time = 0., solution, comparisonData, status,
      reducedMetadata, convergenceData, optimizer, nonLinearResidual,
      kirchhoffResidual, isValidShortcut, telemetry, rejectReason, finalResult},
    model = BuildCriticalQuadraticObjective[d2e];
    If[model === $Failed,
      Return[
        <|
          "Outcome" -> "Unavailable",
          "Telemetry" -> <|
            "QuadraticShortcutAttempted" -> False,
            "QuadraticShortcutAccepted" -> False,
            "QuadraticShortcutNonLinearResidual" -> Missing["NotComputed"],
            "QuadraticShortcutKirchhoffResidual" -> Missing["NotComputed"],
            "QuadraticShortcutSolveTime" -> 0.,
            "QuadraticShortcutRejectReason" -> "ModelUnavailable"
          |>
        |>,
        Module
      ]
    ];
    vars = model["Variables"];
    reducedMetadata = BuildReducedKirchhoffMetadata[model["StateData"]];
    optimizer = If[UseQuadraticCriticalBackendQ[d2e], QuadraticOptimization, ConvexOptimization];
    {time, rawSolution} =
      AbsoluteTiming @ Quiet @ Check[
        optimizer[
          model["Objective"],
          Join[
            Thread[Normal[model["KM"]] . vars == Normal[model["B"]]],
            Thread[vars >= 0]
          ],
          vars
        ],
        $Failed
      ];
    If[!MatchQ[rawSolution, {_Rule ..}],
      telemetry = <|
        "QuadraticShortcutAttempted" -> True,
        "QuadraticShortcutAccepted" -> False,
        "QuadraticShortcutNonLinearResidual" -> Missing["NotComputed"],
        "QuadraticShortcutKirchhoffResidual" -> Missing["NotComputed"],
        "QuadraticShortcutSolveTime" -> time,
        "QuadraticShortcutRejectReason" -> "OptimizerFailed"
      |>;
      Return[<|"Outcome" -> "Unavailable", "Telemetry" -> telemetry|>, Module]
    ];
    solution = AssociationThread[vars, N[vars /. rawSolution]];
    solution = ReconstructUtilities[solution, d2e];
    comparisonData = BuildMonotoneComparisonData[d2e, model["StateData"], solution];
    nonLinearResidual = ComputeMonotoneEquationResidual[d2e, solution];
    comparisonData = Join[comparisonData, <|"NonLinearResidual" -> nonLinearResidual|>];
    status = CheckFlowFeasibility[solution];
    kirchhoffResidual = Lookup[comparisonData, "KirchhoffResidual", Missing["NotAvailable"]];
    isValidShortcut =
      TrueQ[
        status === "Feasible" &&
        NumericQ[nonLinearResidual] &&
        nonLinearResidual <= residualTolerance &&
        NumericQ[kirchhoffResidual] &&
        kirchhoffResidual <= residualTolerance
      ];
    rejectReason =
      Which[
        isValidShortcut, None,
        status =!= "Feasible", "FlowInfeasible",
        !NumericQ[nonLinearResidual], "ResidualNotComputable",
        nonLinearResidual > residualTolerance, "NonLinearResidualExceedsTolerance",
        !NumericQ[kirchhoffResidual], "KirchhoffResidualNotComputable",
        kirchhoffResidual > residualTolerance, "KirchhoffResidualExceedsTolerance",
        True, "Unknown"
      ];
    telemetry = <|
      "QuadraticShortcutAttempted" -> True,
      "QuadraticShortcutAccepted" -> isValidShortcut,
      "QuadraticShortcutNonLinearResidual" -> nonLinearResidual,
      "QuadraticShortcutKirchhoffResidual" -> kirchhoffResidual,
      "QuadraticShortcutSolveTime" -> time,
      "QuadraticShortcutRejectReason" -> rejectReason
    |>;
    If[isValidShortcut,
      convergenceData = <|
        "StopReason" -> "QuadraticCriticalSolve",
        "FinalResidual" -> kirchhoffResidual,
        "Iterations" -> 0,
        "SolveTime" -> time,
        "SeedMethod" -> "QuadraticCriticalSolve"
      |>;
      finalResult = MakeSolverResult[
        "Monotone",
        "Success",
        "Feasible",
        None,
        solution,
        Join[
          comparisonData,
          reducedMetadata,
          <|
            "AssoMonotone" -> solution,
            "Convergence" -> Join[convergenceData, telemetry]
          |>
        ]
      ];
      <|"Outcome" -> "Accepted", "Result" -> finalResult, "Telemetry" -> telemetry|>,
      <|
        "Outcome" -> "Rejected",
        "Telemetry" -> telemetry,
        "CandidateSolution" -> solution
      |>
    ]
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

(* Exact affine reduction x = x0 + N . theta for the Kirchhoff manifold.
   We quotient out directions that are invisible to signed edge flows, so the
   reduced state matches the observable stationary problem much more closely. *)
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

ReducedKirchhoffVector[reduced_Association, theta_?NumberVectorQ] :=
  reduced["BasePoint"] + reduced["BasisMatrix"] . theta;

ReducedKirchhoffAssociation[reduced_Association, theta_?NumberVectorQ] :=
  AssociationThread[reduced["FullVariables"], ReducedKirchhoffVector[reduced, theta]];

(* --- BuildMonotoneField: total cost-to-go field on Kirchhoff variables --- *)
(*
   The stationary system chooses transitions by downstream value, not by local
   edge cost alone. On branching networks, a purely lifted edge-cost field sends
   mass to the locally cheaper branch even when the exit-value equations imply a
   different total cost-to-go. We therefore rebuild the Monotone field from the
   linear u-system induced by the current signed edge costs.
*)

BuildMonotonePairCostAssociation[halfPairs_List, edgeList_List, q_?NumberVectorQ] :=
  Association @ Flatten @ MapThread[
    Function[{pair, edge, flow},
      With[{cost = If[PossibleZeroQ[flow], 0., N @ Cost[Abs[flow], edge]]},
        {pair -> cost, Reverse[pair] -> cost}
      ]
    ],
    {halfPairs, edgeList, q}
  ];

ComputeMonotoneEquationResidual[d2e_Association, solution_Association] :=
  Module[{nlhs, nrhs, diff},
    nlhs = Lookup[d2e, "Nlhs", {}];
    nrhs = Lookup[d2e, "Nrhs", {}];
    If[
      !ListQ[nlhs] || !ListQ[nrhs] || nlhs === {} || nrhs === {} ||
      Length[nlhs] =!= Length[nrhs],
      Return[Missing["NoData"], Module]
    ];
    diff = Quiet @ Check[(nlhs - nrhs) /. solution, $Failed];
    If[diff === $Failed,
      Return[Missing["ComputeError"], Module]
    ];
    diff = Quiet @ Check[N @ diff, $Failed];
    If[diff === $Failed || !VectorQ[diff, NumericQ],
      Missing["NotComputable"],
      Quiet @ Check[N @ Norm[diff, Infinity], Missing["ComputeError"]]
    ]
  ];

ReconstructMonotoneEdgeFlows[solution_Association, d2e_Association] :=
  Module[{entryInRules, flowRules, replacementRules, resolve, edgeFlows},
    entryInRules = Association @ Flatten[ToRules /@ Lookup[d2e, "EqEntryIn", {}]];
    flowRules = Join[
      entryInRules,
      Lookup[d2e, "RuleEntryOut", <||>],
      Lookup[d2e, "RuleExitFlowsIn", <||>],
      Lookup[d2e, "RuleBalanceGatheringFlows", <||>]
    ];
    replacementRules = Join[Normal[solution], Normal[flowRules]];
    resolve[expr_] := FixedPoint[ReplaceAll[#, replacementRules] &, expr, 10];
    edgeFlows = Association @ Cases[
      Map[
        Function[{var},
          Module[{val},
            val = Quiet @ Check[N @ resolve[var], Missing["NotAvailable"]];
            If[NumericQ[val], var -> val, Nothing]
          ]
        ],
        Lookup[d2e, "js", {}]
      ],
      _Rule
    ];
    Join[edgeFlows, solution]
  ];

ReconstructUtilities[solution_Association, d2e_Association] :=
  Module[{solutionWithFlows, stateData, jj, jjVals, halfPairs, edgeList, q,
      pairCosts, valueSystem, utilityAsso},
    solutionWithFlows = ReconstructMonotoneEdgeFlows[solution, d2e];
    stateData = BuildMonotoneStateData[d2e];
    jj = stateData["FullVariables"];
    jjVals = Lookup[solutionWithFlows, jj, Missing["NotAvailable"]];
    If[!ListQ[jjVals] || !VectorQ[jjVals, NumericQ],
      Return[solutionWithFlows, Module]
    ];
    halfPairs = stateData["HalfPairs"];
    edgeList = Lookup[d2e, "edgeList", {}];
    q = If[Length[halfPairs] === 0, {}, N @ (stateData["SignedEdgeMatrix"] . jjVals)];
    pairCosts =
      If[Length[halfPairs] === 0,
        <||>,
        BuildMonotonePairCostAssociation[halfPairs, edgeList, q]
      ];
    valueSystem = BuildMonotoneValueSystem[d2e];
    utilityAsso = Lookup[valueSystem, "StateValueAssociation", (<||> &)][pairCosts];
    If[AssociationQ[utilityAsso],
      Join[utilityAsso, solutionWithFlows],
      solutionWithFlows
    ]
  ];

BuildMonotoneField[d2e_Association] :=
  Module[{stateData, valueSystem, jj, halfPairs, edgeList, S, switching, fieldFn},
    stateData = BuildMonotoneStateData[d2e];
    valueSystem = BuildMonotoneValueSystem[d2e];
    jj = stateData["FullVariables"];
    If[Length[jj] === 0, Return[{jj, (0 * # &)}]];
    halfPairs = stateData["HalfPairs"];
    If[Length[halfPairs] === 0, Return[{jj, (0 * # &)}]];
    edgeList = Lookup[d2e, "edgeList", {}];
    S = stateData["SignedEdgeMatrix"];
    switching = Lookup[d2e, "SwitchingCosts", <||>];
    fieldFn[x_?NumberVectorQ] := Module[{q, pairCosts, stateValues, values},
      q = S . x;
      pairCosts = BuildMonotonePairCostAssociation[halfPairs, edgeList, q];
      stateValues = valueSystem["StateValueAssociation"][pairCosts];
      values = Join[stateValues, <|"PairCosts" -> pairCosts|>];
      MonotoneVariableFieldValue[#, values, switching] & /@ jj
    ];
    {jj, fieldFn}
  ];

BuildReducedMonotoneDynamics[reduced_Association, KM_, B_, cc_, useCache_:True] :=
  Module[{basis, basePoint, dim, edgeOffset, edgeBasis, fullState, edgeState,
      reducedMetric, reducedSolve,
      fullDirection, reducedDirection, residualData},
    basis = reduced["BasisMatrix"];
    basePoint = reduced["BasePoint"];
    dim = reduced["StateDimension"];
    edgeOffset = reduced["SignedEdgeOffset"];
    edgeBasis = reduced["SignedEdgeBasis"];
    fullState[theta_?NumberVectorQ] :=
      If[dim === 0, basePoint, basePoint + basis . theta];
    edgeState[theta_?NumberVectorQ] :=
      If[dim === 0, edgeOffset, edgeOffset + edgeBasis . theta];
    reducedMetric[q_?NumberVectorQ] :=
      If[dim === 0 || Length[q] === 0,
        ConstantArray[0., {dim, dim}],
        Transpose[edgeBasis] . Hess[q] . edgeBasis
      ];
    reducedSolve[q_?NumberVectorQ, rhs_?NumberVectorQ] :=
      Module[{metric, solver},
        metric = reducedMetric[q];
        solver = Quiet @ Check[LinearSolve[metric], $Failed];
        If[solver === $Failed,
          -PseudoInverse[metric] . rhs,
          -solver[rhs]
        ]
      ];
    reducedDirection[theta_?NumberVectorQ] :=
      Module[{x, q},
        x = fullState[theta];
        q = edgeState[theta];
        If[dim === 0, {}, reducedSolve[q, Transpose[basis] . cc[x]]]
      ];
    fullDirection[theta_?NumberVectorQ] :=
      If[dim === 0,
        ConstantArray[0., Length[basePoint]],
        basis . reducedDirection[theta]
      ];
    residualData[theta_?NumberVectorQ] :=
      Module[{x, q, direction, reducedDir, stationarity, fullStationarity, kirchhoff},
        x = fullState[theta];
        q = edgeState[theta];
        reducedDir = reducedDirection[theta];
        direction = fullDirection[theta];
        stationarity = If[Length[reducedDir] === 0, 0., N @ Norm[reducedDir, Infinity]];
        fullStationarity = If[Length[direction] === 0, 0., N @ Norm[direction, Infinity]];
        kirchhoff = N @ Norm[KM . x - B, Infinity];
        <|
          "Solution" -> AssociationThread[reduced["FullVariables"], x],
          "StationarityResidual" -> stationarity,
          "FullStationarityResidual" -> fullStationarity,
          "KirchhoffResidual" -> kirchhoff,
          "FinalResidual" -> Max[stationarity, kirchhoff],
          "MinFlow" -> If[Length[q] === 0, Infinity, N @ Min[q]]
        |>
      ];
    <|
      "StateDimension" -> dim,
      "FullState" -> fullState,
      "ReducedDirection" -> reducedDirection,
      "ResidualData" -> residualData
    |>
  ];

(* --- MonotoneSolver --- *)

MonotoneSolverFromData[Data_, opts:OptionsPattern[]] :=
	Module[{d2e, potentialFunction, congestionExponentFunction, interactionFunction},
		potentialFunction = OptionValue["PotentialFunction"];
		congestionExponentFunction = OptionValue["CongestionExponentFunction"];
		interactionFunction = OptionValue["InteractionFunction"];
		d2e = WithHamiltonianFunctions[
			potentialFunction, congestionExponentFunction, interactionFunction,
			DataToEquations[Data]
		];
		MonotoneSolver[d2e, opts]
	];

MonotoneSolver[d2e_, opts:OptionsPattern[]] :=
	Module[{potentialFunction, congestionExponentFunction, interactionFunction},
		potentialFunction = OptionValue["PotentialFunction"];
		congestionExponentFunction = OptionValue["CongestionExponentFunction"];
		interactionFunction = OptionValue["InteractionFunction"];
		WithHamiltonianFunctions[
			potentialFunction,
			congestionExponentFunction,
			interactionFunction,
				Module[{stateData, reducedState, B, KM, jj, cc, x0, result, solution, odeResult,
				    missingSolution, comparisonData, reducedMetadata, convergenceData,
                    resultKind, message, residualTolerance, status, nonLinearResidual,
                    residualViolationQ, shortcutCandidate, shortcutOutcome,
                    shortcutTelemetry, defaultShortcutTelemetry, candidateSolution,
                    warmStartVector, seedMethod, warmStartFloor},
					missingSolution = Missing["NotAvailable"];
                    residualTolerance = N @ OptionValue["ResidualTolerance"];
                    warmStartFloor = 10^-6;
                    seedMethod = "InteriorFindInstance";
					stateData = BuildMonotoneStateData[d2e];
					B = stateData["B"];
					KM = stateData["KM"];
					jj = stateData["FullVariables"];
					(* Guard: skip MonotoneSolver for degenerate cases with no flow variables *)
					If[Length[jj] === 0,
						Message[MonotoneSolver::degenerate];
                        comparisonData = BuildSolverComparisonData[d2e, missingSolution];
                        convergenceData = <|
                            "StopReason" -> "DegenerateCase",
                            "FinalResidual" -> Missing["NotApplicable"],
                            "Iterations" -> 0,
                            "SolveTime" -> 0.,
                            "SeedMethod" -> "NotApplicable"
                        |>;
						Return[
							MakeSolverResult[
								"Monotone",
								"Degenerate",
								Missing["NotApplicable"],
								"DegenerateCase",
								missingSolution,
								Join[comparisonData, <|
                                    "AssoMonotone" -> missingSolution,
                                    "Convergence" -> convergenceData
                                |>]
								]
					]
					];
					(* Build the lifted edge-cost field *)
                    defaultShortcutTelemetry = <|
                        "QuadraticShortcutAttempted" -> False,
                        "QuadraticShortcutAccepted" -> False,
                        "QuadraticShortcutNonLinearResidual" -> Missing["NotComputed"],
                        "QuadraticShortcutKirchhoffResidual" -> Missing["NotComputed"],
                        "QuadraticShortcutSolveTime" -> 0.,
                        "QuadraticShortcutRejectReason" -> "ModelUnavailable"
                    |>;
                    shortcutCandidate = TryCriticalQuadraticMonotoneSolve[d2e, residualTolerance];
                    shortcutOutcome = Lookup[shortcutCandidate, "Outcome", "Unavailable"];
                    shortcutTelemetry = Lookup[shortcutCandidate, "Telemetry", defaultShortcutTelemetry];
                    candidateSolution = Lookup[shortcutCandidate, "CandidateSolution", Missing["NotAvailable"]];
                    If[
                        shortcutOutcome === "Accepted" &&
                        AssociationQ[shortcutCandidate] &&
                        KeyExistsQ[shortcutCandidate, "Result"],
                        Return[shortcutCandidate["Result"]]
                    ];
					{jj, cc} = BuildMonotoneField[d2e];
                    (* Warm start from rejected quadratic candidate when available *)
                    warmStartVector =
                        If[AssociationQ[candidateSolution],
                            Lookup[candidateSolution, jj, Missing["NotAvailable"]],
                            Missing["NotAvailable"]
                        ];
                    If[ListQ[warmStartVector] && VectorQ[warmStartVector, NumericQ],
                        x0 = N @ Map[Max[#, warmStartFloor] &, warmStartVector];
                        seedMethod = "QuadraticWarmStart";
                        ,
                        (* Find numerically interior feasible seed *)
                        result = Quiet @ Check[
                            First @ FindInstance[KM . jj == B && And @@ ((# >= 10^-8) & /@ jj), jj, Reals],
                        $Failed
                    ];
                        If[result === $Failed,
                            Message[MonotoneSolver::seedfail];
                            comparisonData = BuildSolverComparisonData[d2e, missingSolution];
                            convergenceData = Join[
                                <|
                                    "StopReason" -> "SeedFindInstanceFailed",
                                    "FinalResidual" -> Missing["NotApplicable"],
                                    "Iterations" -> 0,
                                    "SolveTime" -> 0.,
                                    "SeedMethod" -> seedMethod
                                |>,
                                shortcutTelemetry
                            ];
                            Return[
                                MakeSolverResult[
                                    "Monotone",
                                    "Failure",
                                    Missing["NotApplicable"],
                                    "SeedFindInstanceFailed",
                                    missingSolution,
                                    Join[comparisonData, <|
                                        "AssoMonotone" -> missingSolution,
                                        "Convergence" -> convergenceData
                                    |>]
                                    ]
                        ]
                    ];
                        x0 = N[jj /. result];
                    ];
                    shortcutTelemetry = Join[shortcutTelemetry, <|"SeedMethod" -> seedMethod|>];
                    reducedState = BuildReducedKirchhoffCoordinates[d2e, x0];
				(* No edges → zero cost field → feasible point is already the solution *)
				odeResult = If[Length[Lookup[d2e, "edgeList", {}]] === 0,
					<|
                            "Solution" -> AssociationThread[jj, x0],
                            "Convergence" -> <|
                                "StopReason" -> "ResidualToleranceMet",
                                "FinalResidual" -> 0.,
                                "Iterations" -> 0,
                                "SolveTime" -> 0.,
                                "SeedMethod" -> seedMethod
                            |>
                        |>,
						MonotoneSolverODE[reducedState, KM, B, cc, FilterRules[{opts}, Options[MonotoneSolverODE]]]
					];
                    solution = Lookup[odeResult, "Solution", missingSolution];
                    convergenceData = Join[
                        Lookup[odeResult, "Convergence",
                            <|
                                "StopReason" -> "ODEFailure",
                                "FinalResidual" -> Missing["NotAvailable"],
                                "Iterations" -> 0,
                                "SolveTime" -> 0.,
                                "SeedMethod" -> seedMethod
                            |>
                        ],
                        shortcutTelemetry
                    ];
                    If[MissingQ[solution],
                        Message[MonotoneSolver::odefail];
                        comparisonData = BuildSolverComparisonData[d2e, missingSolution];
                        Return[
                            MakeSolverResult[
                                "Monotone",
                                "Failure",
                                Missing["NotApplicable"],
                                "ODEFailure",
                                missingSolution,
                                Join[comparisonData, <|
                                    "AssoMonotone" -> missingSolution,
                                    "Convergence" -> convergenceData
                                |>]
                            ]
                        ]
                    ];
					solution = ReconstructUtilities[solution, d2e];
					comparisonData = BuildSolverComparisonData[d2e, solution];
                    nonLinearResidual = ComputeMonotoneEquationResidual[d2e, solution];
                    comparisonData = Join[comparisonData, <|"NonLinearResidual" -> nonLinearResidual|>];
                    status = CheckFlowFeasibility[solution];
                    residualViolationQ =
                        NumericQ[nonLinearResidual] && nonLinearResidual > residualTolerance;
                    resultKind =
                        If[
                            !residualViolationQ &&
                            status === "Feasible" &&
                            NumericQ[Lookup[convergenceData, "FinalResidual", Missing["NotAvailable"]]] &&
                            Lookup[convergenceData, "FinalResidual", Infinity] <= residualTolerance &&
                            NumericQ[Lookup[comparisonData, "KirchhoffResidual", Missing["NotAvailable"]]] &&
                            Lookup[comparisonData, "KirchhoffResidual", Infinity] <= residualTolerance,
                            "Success",
                            "NonConverged"
                        ];
                    If[residualViolationQ, status = "Infeasible"];
                    message =
                        Which[
                            resultKind === "Success", None,
                            residualViolationQ, "ResidualExceedsTolerance",
                            status =!= "Feasible", "InfeasibleProjectedSolution",
                            True, Lookup[convergenceData, "StopReason", "NonConverged"]
                        ];
					reducedMetadata = <|
						"ReducedStateDimension" -> reducedState["StateDimension"],
						"FullStateDimension" -> reducedState["FullDimension"],
						"CostInvisibleDimension" -> reducedState["CostInvisibleDimension"]
					|>;
					MakeSolverResult[
						"Monotone",
						resultKind,
						status,
						message,
						solution,
						Join[comparisonData, reducedMetadata, <|
                            "AssoMonotone" -> solution,
                            "Convergence" -> convergenceData
                        |>]
					]
			]
		]
	];

MonotoneSolverODE[reduced_Association, KM_, B_, cc_, opts:OptionsPattern[]] :=
	Module[{residualTolerance, maxTime, maxSteps, useCache, dynamics, dim, theta0,
        iterations = 0, stopReason = None, solveTime, sol, messages, endTime, thetaEnd,
        finalData, convergenceData},
		residualTolerance = N @ OptionValue["ResidualTolerance"];
        maxTime = N @ OptionValue["MaxTime"];
		maxSteps = OptionValue["MaxSteps"];
		useCache = OptionValue["UseCachedProjection"];
        dynamics = BuildReducedMonotoneDynamics[reduced, KM, B, cc, useCache];
        dim = dynamics["StateDimension"];
        theta0 = ConstantArray[0., dim];
        If[dim === 0,
            finalData = dynamics["ResidualData"][theta0];
            convergenceData = <|
                "StopReason" -> If[finalData["FinalResidual"] <= residualTolerance, "ResidualToleranceMet", "NoDegreesOfFreedom"],
                "FinalResidual" -> finalData["FinalResidual"],
                "Iterations" -> 0,
                "SolveTime" -> 0.,
                "SeedMethod" -> "InteriorFindInstance"
            |>;
            Return[<|"Solution" -> finalData["Solution"], "Convergence" -> convergenceData|>]
        ];
        {solveTime, sol} = AbsoluteTiming[
            Block[{$MessageList = {}},
                iterations = 0;
                sol = Quiet[
                    NDSolveValue[
                        {
                            theta'[t] == dynamics["ReducedDirection"][theta[t]],
                            theta[0] == theta0,
                            WhenEvent[
                                dynamics["ResidualData"][theta[t]]["FinalResidual"] <= residualTolerance,
                                {stopReason = "ResidualToleranceMet", "StopIntegration"}
                            ],
                            WhenEvent[
                                dynamics["ResidualData"][theta[t]]["MinFlow"] <= 10^-10,
                                {If[stopReason === None, stopReason = "BoundaryHit"], "StopIntegration"}
                            ]
                        },
                        theta,
                        {t, 0, maxTime},
                        Method -> "BDF",
                        MaxSteps -> maxSteps,
                        StepMonitor :> (iterations++)
                    ],
                    {NDSolveValue::mxst, NDSolveValue::nbnum1}
                ];
                messages = $MessageList;
                sol
            ]
        ];
        If[Head[sol] =!= InterpolatingFunction,
            Return[
                <|
                    "Solution" -> Missing["NotAvailable"],
                    "Convergence" -> <|
                        "StopReason" -> "ODEFailure",
                        "FinalResidual" -> Missing["NotAvailable"],
                        "Iterations" -> iterations,
                        "SolveTime" -> solveTime,
                        "SeedMethod" -> "InteriorFindInstance"
                    |>
                |>
            ]
        ];
        endTime = Last @ Last[sol["Domain"]];
        thetaEnd = N @ sol[endTime];
        finalData = dynamics["ResidualData"][thetaEnd];
        If[stopReason === None,
            stopReason = Which[
                NumericQ[finalData["FinalResidual"]] && finalData["FinalResidual"] <= residualTolerance, "ResidualToleranceMet",
                IntegerQ[maxSteps] && iterations >= maxSteps, "MaxStepsReached",
                MemberQ[messages, HoldForm[MessageName[NDSolveValue, mxst]]], "MaxStepsReached",
                endTime >= maxTime - 10^-8, "MaxTimeReached",
                True, "Stopped"
            ]
        ];
        convergenceData = <|
            "StopReason" -> stopReason,
            "FinalResidual" -> finalData["FinalResidual"],
            "Iterations" -> iterations,
            "SolveTime" -> solveTime,
            "SeedMethod" -> "InteriorFindInstance"
        |>;
        <|"Solution" -> finalData["Solution"], "Convergence" -> convergenceData|>
	];

End[];
