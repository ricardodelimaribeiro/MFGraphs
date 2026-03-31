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
that caches the PseudoInverse result and reuses it when x has not changed significantly.
cache must be a symbol holding a 1-element list {Null} or {<|\"x\" -> ..., \"pi\" -> ...|>}.
The function has the HoldAll attribute so that cache is passed by reference.";

MonotoneSolverFromData::usage =
"MonotoneSolverFromData[Data] solves the MFG problem from raw Data using the monotone operator method.
Options: \"TimeSteps\" (default 100), \"UseCachedGradient\" (default True), \
\"ReturnShape\" (default \"Legacy\"; use \"Standard\" for a normalized solver-result association), \
\"PotentialFunction\", \"CongestionExponentFunction\", and \"InteractionFunction\" \
(default Automatic, meaning use the current global MFGraphs` definitions of V, alpha, and g).";

MonotoneSolver::usage =
"MonotoneSolver[d2e] solves the MFG problem using the monotone operator method.
Options: \"TimeSteps\" (default 100), \"UseCachedGradient\" (default True), \
\"ReturnShape\" (default \"Legacy\"; use \"Standard\" for a normalized solver-result association), \
\"PotentialFunction\", \"CongestionExponentFunction\", and \"InteractionFunction\" \
(default Automatic, meaning use the current global MFGraphs` definitions of V, alpha, and g).";

MonotoneSolverODE::usage =
"MonotoneSolverODE[x0, KM, jj, cc] solves the gradient flow ODE from initial condition x0.
Options: \"TimeSteps\" (default 100), \"UseCachedGradient\" (default True).";

MonotoneSolver::degenerate =
"Degenerate case: no flow variables were found, so the monotone solver was skipped.";

MonotoneSolver::seedfail =
"Failed to find an interior feasible seed for the monotone ODE solve.";

Options[MonotoneSolverFromData] = {
    "TimeSteps" -> 100,
    "UseCachedGradient" -> True,
    "ReturnShape" -> "Legacy",
    "PotentialFunction" -> Automatic,
    "CongestionExponentFunction" -> Automatic,
    "InteractionFunction" -> Automatic
};
Options[MonotoneSolver] = Options[MonotoneSolverFromData];
Options[MonotoneSolverODE] = {"TimeSteps" -> 100, "UseCachedGradient" -> True};

Begin["`Private`"];

(* HoldAll so that the cache argument (a mutable list) is passed by reference *)
SetAttributes[CachedGradientProjection, HoldAll];

(* --- Matrix utilities --- *)

Hess[j_?NumberVectorQ] := DiagonalMatrix[1/# & /@ j];

InverseHessian[j_?NumberVectorQ] := DiagonalMatrix[j];

NumberMatrixQ[A_] := NumberVectorQ@Flatten@A;

HessianSandwich[x_?NumberVectorQ, A_, At_] := A . InverseHessian[x] . At;
HessianSandwich[x_?NumberVectorQ, A_] := A . InverseHessian[x] . Transpose[A];

(* --- GradientProjection: projected gradient operator --- *)

GradientProjection[x_?NumberVectorQ, A_, dim_, At_] :=
  InverseHessian[x] . (IdentityMatrix[dim] - At . PseudoInverse[HessianSandwich[x, A, At]] . A . InverseHessian[x]);

GradientProjection[x_?NumberVectorQ, A_] :=
	Module[{dim, At},
		dim = Length[x];
		At = Transpose[A];
		InverseHessian[x] . (IdentityMatrix[dim] - At . PseudoInverse[HessianSandwich[x, A, At]] . A . InverseHessian[x])
	];

(* --- CachedGradientProjection: caches PseudoInverse when x changes slowly --- *)
(* HoldAll keeps cache as the caller's symbol so part-assignment works. *)
(* Other arguments are force-evaluated via Module initializers. *)

CachedGradientProjection[x_, KM_, dim_, At_, cache_, tol_:10^-6] :=
  Module[{xVal = x, KMVal = KM, dimVal = dim, AtVal = At, tolVal = tol,
          invH, AinvHAt, pi},
    invH = InverseHessian[xVal];
    If[cache[[1]] =!= Null && Norm[xVal - cache[[1]]["x"], Infinity] < tolVal,
      (* Reuse cached PseudoInverse *)
      invH . (IdentityMatrix[dimVal] - AtVal . cache[[1]]["pi"] . KMVal . invH),
      (* Recompute PseudoInverse and cache *)
      AinvHAt = KMVal . invH . AtVal;
      pi = PseudoInverse[AinvHAt];
      cache[[1]] = <|"x" -> xVal, "pi" -> pi|>;
      invH . (IdentityMatrix[dimVal] - AtVal . pi . KMVal . invH)
    ]
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

(* Exact affine reduction x = x0 + N . theta for the Kirchhoff manifold.
   This does not yet implement the paper's edge-current reduction, but it gives
   a minimal state dimension and lossless reconstruction of the full j[...] state. *)
BuildReducedKirchhoffCoordinates[d2e_Association, basePoint_:Automatic] :=
  Module[{stateData, KM, jj, S, nullBasis, basisMatrix, dim, baseVec, edgeOffset,
      edgeBasis, invisibleDim, signedEdgeRows},
    stateData = BuildMonotoneStateData[d2e];
    KM = stateData["KM"];
    jj = stateData["FullVariables"];
    S = stateData["SignedEdgeMatrix"];
    nullBasis = NullSpace[Normal[KM]];
    basisMatrix = If[nullBasis === {},
      ConstantArray[0., {Length[jj], 0}],
      N @ Transpose[nullBasis]
    ];
    dim = If[MatrixQ[basisMatrix], Last @ Dimensions[basisMatrix], 0];
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
    invisibleDim = If[signedEdgeRows === 0,
      Length[NullSpace[Normal[KM]]],
      Length[NullSpace[Join[Normal[KM], Normal[S]]]]
    ];
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

(* --- BuildMonotoneField: lifted edge-cost field on Kirchhoff variables --- *)
(* 
   This function constructs the vector field for the HRF method by "lifting" 
   microscopic edge costs to the reduced transition-flow space (the ϑ vector 
   from the HRF paper).
   
   Returns {jj, fieldFn} where fieldFn[x_?NumberVectorQ] gives S^T . c(S . x). 
   Here S is the sparsity matrix mapping transition flows to signed edge flows.
*)

BuildMonotoneField[d2e_Association] :=
  Module[{stateData, jj, halfPairs, S, fieldFn},
    stateData = BuildMonotoneStateData[d2e];
    jj = stateData["FullVariables"];
    If[Length[jj] === 0, Return[{jj, (0 * # &)}, Module]];
    halfPairs = stateData["HalfPairs"];
    If[Length[halfPairs] === 0, Return[{jj, (0 * # &)}, Module]];
    S = stateData["SignedEdgeMatrix"];
    fieldFn[x_?NumberVectorQ] := Module[{q, c},
      q = S . x;
      c = MapThread[
        If[PossibleZeroQ[#1], 0., Sign[#1] Cost[#1, #2]] &,
        {q, halfPairs}
      ];
      Transpose[S] . c
    ];
    {jj, fieldFn}
  ];

(* --- MonotoneSolver --- *)

MonotoneSolverFromData[Data_, opts:OptionsPattern[]] :=
	Module[{d2e},
		d2e = DataToEquations[Data];
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
				Module[{stateData, reducedState, B, KM, jj, cc, x0, result, returnShape, solution,
				    missingSolution, comparisonData, reducedMetadata},
					returnShape = OptionValue["ReturnShape"];
					missingSolution = Missing["NotAvailable"];
					stateData = BuildMonotoneStateData[d2e];
					B = stateData["B"];
					KM = stateData["KM"];
					jj = stateData["FullVariables"];
					(* Guard: skip MonotoneSolver for degenerate cases with no flow variables *)
					If[Length[jj] === 0,
						Message[MonotoneSolver::degenerate];
						Return[
							If[returnShape === "Standard",
								comparisonData = BuildSolverComparisonData[d2e, missingSolution];
								MakeSolverResult[
									"Monotone",
									"Degenerate",
									Missing["NotApplicable"],
									"DegenerateCase",
									missingSolution,
									Join[comparisonData, <|"AssoMonotone" -> missingSolution|>]
								],
								<|"Message" -> "Degenerate case"|>
							],
						Module
					]
				];
				(* Build the lifted edge-cost field *)
				{jj, cc} = BuildMonotoneField[d2e];
				(* Find numerically interior feasible seed *)
				result = Quiet @ Check[
					First @ FindInstance[KM . jj == B && And @@ ((# >= 10^-8) & /@ jj), jj, Reals],
					$Failed
				];
					If[result === $Failed,
						Message[MonotoneSolver::seedfail];
						Return[
							If[returnShape === "Standard",
								comparisonData = BuildSolverComparisonData[d2e, missingSolution];
								MakeSolverResult[
									"Monotone",
									"Failure",
									Missing["NotApplicable"],
									"SeedFindInstanceFailed",
									missingSolution,
									Join[comparisonData, <|"AssoMonotone" -> missingSolution|>]
								],
								Null
							],
						Module
					]
				];
				x0 = N[jj /. result];
				(* No edges → zero cost field → feasible point is already the solution *)
				solution = If[Length[Lookup[d2e, "edgeList", {}]] === 0,
					AssociationThread[jj, x0],
						MonotoneSolverODE[x0, KM, jj, cc, FilterRules[{opts}, Options[MonotoneSolverODE]]]
					];
					If[returnShape === "Standard",
						reducedState = BuildReducedKirchhoffCoordinates[d2e, x0];
						comparisonData = BuildSolverComparisonData[d2e, solution];
						reducedMetadata = <|
							"ReducedStateDimension" -> reducedState["StateDimension"],
							"FullStateDimension" -> reducedState["FullDimension"],
							"CostInvisibleDimension" -> reducedState["CostInvisibleDimension"]
						|>;
						MakeSolverResult[
							"Monotone",
							"Success",
							Missing["NotApplicable"],
							None,
							solution,
							Join[comparisonData, reducedMetadata, <|"AssoMonotone" -> solution|>]
						],
						solution
					]
			]
		]
	];

MonotoneSolverODE[x0_, KM_, jj_, cc_, opts:OptionsPattern[]] :=
	Module[{At, dim, sol, x, xx, n, useCache, piCache, gradFn},
		n = OptionValue["TimeSteps"];
		useCache = OptionValue["UseCachedGradient"];
		At = Transpose[KM];
		dim = Length[x0];

		If[useCache,
			(* CachedGradientProjection has HoldAll; piCache is passed by reference *)
			piCache = {Null};
			gradFn[xv_?NumberVectorQ] := CachedGradientProjection[xv, KM, dim, At, piCache],
			(* Use original uncached version *)
			gradFn[xv_?NumberVectorQ] := GradientProjection[xv, KM, dim, At]
		];

		sol = NDSolve[
			x'[t] == -gradFn[x[t]] . cc[x[t]] && x[0] == x0,
			x[t], {t, 0, n}, Method -> "BDF"];
		xx = Table[x[t] /. sol[[1]], {t, 1, n, 1}];
		AssociationThread[jj, Last[xx]]
	];

End[];
