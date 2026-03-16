(* Wolfram Language package *)
(* Monotone operator solver for Mean Field Games on networks *)

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
"CachedGradientProjection[x, K, dim, At, cache] is a version of GradientProjection
that caches the PseudoInverse result and reuses it when x has not changed significantly.
cache must be a symbol holding a 1-element list {Null} or {<|\"x\" -> ..., \"pi\" -> ...|>}.
The function has the HoldAll attribute so that cache is passed by reference.";

MonotoneSolverFromData::usage =
"MonotoneSolverFromData[Data] solves the MFG problem from raw Data using the monotone operator method.
Options: \"TimeSteps\" (default 100), \"UseCachedGradient\" (default True).";

MonotoneSolver::usage =
"MonotoneSolver[d2e] solves the MFG problem using the monotone operator method.
Options: \"TimeSteps\" (default 100), \"UseCachedGradient\" (default True).";

MonotoneSolverODE::usage =
"MonotoneSolverODE[x0, K, jj, cc] solves the gradient flow ODE from initial condition x0.
Options: \"TimeSteps\" (default 100), \"UseCachedGradient\" (default True).";

Options[MonotoneSolverFromData] = {"TimeSteps" -> 100, "UseCachedGradient" -> True};
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

CachedGradientProjection[x_, K_, dim_, At_, cache_, tol_:10^-6] :=
  Module[{xVal = x, KVal = K, dimVal = dim, AtVal = At, tolVal = tol,
          invH, AinvHAt, pi},
    invH = InverseHessian[xVal];
    If[cache[[1]] =!= Null && Norm[xVal - cache[[1]]["x"], Infinity] < tolVal,
      (* Reuse cached PseudoInverse *)
      invH . (IdentityMatrix[dimVal] - AtVal . cache[[1]]["pi"] . KVal . invH),
      (* Recompute PseudoInverse and cache *)
      AinvHAt = KVal . invH . AtVal;
      pi = PseudoInverse[AinvHAt];
      cache[[1]] = <|"x" -> xVal, "pi" -> pi|>;
      invH . (IdentityMatrix[dimVal] - AtVal . pi . KVal . invH)
    ]
  ];

(* --- BuildMonotoneField: lifted edge-cost field on Kirchhoff variables --- *)
(* Returns {jj, fieldFn} where fieldFn[x_?NumberVectorQ] gives S^T . c(S . x). *)

BuildMonotoneField[d2e_Association] :=
  Module[{B, K, cost, jj, halfPairs, signedFlows, rules, substituted, S, fieldFn},
    {B, K, cost, jj} = GetKirchhoffMatrix[d2e];
    If[Length[jj] === 0, Return[{jj, (0 * # &)}, Module]];
    halfPairs = List @@@ Lookup[d2e, "edgeList", {}];
    If[Length[halfPairs] === 0, Return[{jj, (0 * # &)}, Module]];
    signedFlows = Lookup[d2e, "SignedFlows", <||>];
    rules = Join[
      Lookup[d2e, "RuleBalanceGatheringFlows", <||>],
      Lookup[d2e, "RuleExitFlowsIn", <||>],
      Lookup[d2e, "RuleEntryOut", <||>]
    ];
    (* Substitute rules into signed flows to get expressions in jj *)
    substituted = (signedFlows[#] & /@ halfPairs) /. rules;
    (* S is the sparse coefficient matrix: signed_edge_flows = S . jj *)
    S = Last @ CoefficientArrays[substituted, jj];
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
	Module[{B, K, cost, jj, cc, x0, result},
		{B, K, cost, jj} = GetKirchhoffMatrix[d2e];
		(* Guard: skip MonotoneSolver for degenerate cases with no flow variables *)
		If[Length[jj] === 0,
			MFGPrint["MonotoneSolver: Degenerate case (no flow variables), skipping."];
			Return[<|"Message" -> "Degenerate case"|>, Module]
		];
		(* Build the lifted edge-cost field *)
		{jj, cc} = BuildMonotoneField[d2e];
		(* Find numerically interior feasible seed *)
		result = Quiet @ Check[
			First @ FindInstance[K . jj == B && And @@ ((# >= 10^-8) & /@ jj), jj, Reals],
			$Failed
		];
		If[result === $Failed,
			MFGPrint["MonotoneSolver: FindInstance failed, returning null solution."];
			Return[Null, Module]
		];
		x0 = N[jj /. result];
		(* No edges → zero cost field → feasible point is already the solution *)
		If[Length[Lookup[d2e, "edgeList", {}]] === 0,
			Return[AssociationThread[jj, x0], Module]
		];
		MonotoneSolverODE[x0, K, jj, cc, FilterRules[{opts}, Options[MonotoneSolverODE]]]
	];

MonotoneSolverODE[x0_, K_, jj_, cc_, opts:OptionsPattern[]] :=
	Module[{At, dim, sol, x, xx, n, useCache, piCache, gradFn},
		n = OptionValue["TimeSteps"];
		useCache = OptionValue["UseCachedGradient"];
		At = Transpose[K];
		dim = Length[x0];

		If[useCache,
			(* CachedGradientProjection has HoldAll; piCache is passed by reference *)
			piCache = {Null};
			gradFn[xv_?NumberVectorQ] := CachedGradientProjection[xv, K, dim, At, piCache],
			(* Use original uncached version *)
			gradFn[xv_?NumberVectorQ] := GradientProjection[xv, K, dim, At]
		];

		sol = NDSolve[
			x'[t] == -gradFn[x[t]] . cc[x[t]] && x[0] == x0,
			x[t], {t, 0, n}, Method -> "BDF"];
		xx = Table[x[t] /. sol[[1]], {t, 1, n, 1}];
		AssociationThread[jj, Last[xx]]
	];

End[];
