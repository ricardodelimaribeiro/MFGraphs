(* Wolfram Language package *)
(* Monotone operator solver for Mean Field Games on networks *)

(* --- Matrix utilities --- *)

Hess::usage =
"Hess[j] returns a numeric diagonal matrix with the reciprocals of the elements in the (numeric) vector j.";
Hess[j_?NumberVectorQ] := DiagonalMatrix[1/# & /@ j];

InverseHessian::usage =
"InverseHessian[j] returns a diagonal matrix with the (numeric) vector j. This is the inverse of the matrix Hess[j].";
InverseHessian[j_?NumberVectorQ] := DiagonalMatrix[j];

NumberMatrixQ::usage =
"NumberMatrixQ[A] returns True if the elements of the matrix A are numeric.";
NumberMatrixQ[A_] := NumberVectorQ@Flatten@A;

HessianSandwich::usage =
"HessianSandwich[j, A, At] returns the product A . InverseHessian[j] . At.";
HessianSandwich[x_?NumberVectorQ, A_, At_] := A . InverseHessian[x] . At;
HessianSandwich[x_?NumberVectorQ, A_] := A . InverseHessian[x] . Transpose[A];

(* --- GradientProjection: projected gradient operator --- *)

GradientProjection::usage =
"GradientProjection[x, A, dim, At] returns the projected gradient operator matrix.
This is InverseHessian[x] . (I - At . PseudoInverse[A . InverseHessian[x] . At] . A . InverseHessian[x]).";

GradientProjection[x_?NumberVectorQ, A_, dim_, At_] :=
  InverseHessian[x] . (IdentityMatrix[dim] - At . PseudoInverse[HessianSandwich[x, A, At]] . A . InverseHessian[x]);

GradientProjection[x_?NumberVectorQ, A_] :=
	Module[{dim, At},
		dim = Length[x];
		At = Transpose[A];
		InverseHessian[x] . (IdentityMatrix[dim] - At . PseudoInverse[HessianSandwich[x, A, At]] . A . InverseHessian[x])
	];

(* --- CachedGradientProjection: caches PseudoInverse when x changes slowly --- *)

CachedGradientProjection::usage =
"CachedGradientProjection[x, K, dim, At, cache] is a version of GradientProjection
that caches the PseudoInverse result and reuses it when x has not changed significantly.
cache should be a held Symbol containing <|\"x\" -> ..., \"pi\" -> ...|> or Null.";

CachedGradientProjection[x_?NumberVectorQ, K_, dim_, At_, cache_Symbol, tol_:10^-6] :=
  Module[{invH, AinvHAt, pi},
    invH = InverseHessian[x];
    If[cache =!= Null && Norm[x - cache["x"], Infinity] < tol,
      (* Reuse cached PseudoInverse *)
      invH . (IdentityMatrix[dim] - At . cache["pi"] . K . invH),
      (* Recompute PseudoInverse and cache *)
      AinvHAt = K . invH . At;
      pi = PseudoInverse[AinvHAt];
      cache = <|"x" -> x, "pi" -> pi|>;
      invH . (IdentityMatrix[dim] - At . pi . K . invH)
    ]
  ];

(* --- MonotoneSolver --- *)

Options[MonotoneSolverFromData] = {"TimeSteps" -> 100, "UseCachedGradient" -> True};
Options[MonotoneSolver] = Options[MonotoneSolverFromData];

MonotoneSolverFromData::usage =
"MonotoneSolverFromData[Data] solves the MFG problem from raw Data using the monotone operator method.
Options: \"TimeSteps\" (default 100), \"UseCachedGradient\" (default True).";
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
		cc[j_?NumberVectorQ] := cost[j];
		(* Guard: check if FindInstance succeeds *)
		result = Quiet @ Check[
			First@FindInstance[K . jj == B && And @@ ((# > 0) & /@ jj), jj],
			$Failed
		];
		If[result === $Failed,
			MFGPrint["MonotoneSolver: FindInstance failed, returning null solution."];
			Return[Null, Module]
		];
		x0 = jj /. result;
		MonotoneSolverODE[x0, K, jj, cc, FilterRules[{opts}, Options[MonotoneSolverODE]]]
	];

Options[MonotoneSolverODE] = {"TimeSteps" -> 100, "UseCachedGradient" -> True};

MonotoneSolverODE::usage =
"MonotoneSolverODE[x0, K, jj, cc] solves the gradient flow ODE from initial condition x0.
Options: \"TimeSteps\" (default 100), \"UseCachedGradient\" (default True).";
MonotoneSolverODE[x0_, K_, jj_, cc_, opts:OptionsPattern[]] :=
	Module[{At, dim, sol, x, xx, n, useCache, piCache, gradFn},
		n = OptionValue["TimeSteps"];
		useCache = OptionValue["UseCachedGradient"];
		At = Transpose[K];
		dim = Length[x0];

		If[useCache,
			(* Use cached version to avoid redundant PseudoInverse calls *)
			piCache = Null;
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
