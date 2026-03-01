(* Wolfram Language package *)
(* Monotone operator solver for Mean Field Games on networks *)

(* --- Matrix utilities --- *)

Hess::usage =
"Hess[j] returns a numeric diagonal matrix with the reciprocals of the elements in the (numeric) vector j.";
Hess[j_?NumberVectorQ] := DiagonalMatrix[1/# & /@ j];

H2[j_?NumberVectorQ] := DiagonalMatrix[1/#^2 & /@ j];

H3[j_?NumberVectorQ] := DiagonalMatrix[1/#^(1/4) & /@ j];

InvH::usage =
"InvH[j] returns a diagonal matrix with the (numeric) vector j. This is the inverse of the matrix Hess[j].";
InvH[j_?NumberVectorQ] := DiagonalMatrix[j];

NumberMatrixQ::usage =
"NumberMatrixQ[A] returns True if the elements of the matrix A are numeric.";
NumberMatrixQ[A_] := NumberVectorQ@Flatten@A;

InvAHA::usage =
"InvAHA[j, A, At] returns the product A . InvH[j] . At.";
InvAHA[x_?NumberVectorQ, A_, At_] := A . InvH[x] . At;
InvAHA[x_?NumberVectorQ, A_] := A . InvH[x] . Transpose[A];

(* --- GradientProjection: projected gradient operator --- *)

GradientProjection::usage =
"GradientProjection[x, A, dim, At] returns the projected gradient operator matrix.
This is InvH[x] . (I - At . PseudoInverse[A . InvH[x] . At] . A . InvH[x]).";

GradientProjection[x_?NumberVectorQ, A_, dim_, At_] :=
  InvH[x] . (IdentityMatrix[dim] - At . PseudoInverse[InvAHA[x, A, At]] . A . InvH[x]);

GradientProjection[x_?NumberVectorQ, A_] :=
	Module[{dim, At},
		dim = Length[x];
		At = Transpose[A];
		InvH[x] . (IdentityMatrix[dim] - At . PseudoInverse[InvAHA[x, A, At]] . A . InvH[x])
	];

(* Backward-compatible alias *)
GG = GradientProjection;

(* --- MonotoneSolver --- *)

Options[MonotoneSolverFromData] = {"TimeSteps" -> 100};
Options[MonotoneSolver] = Options[MonotoneSolverFromData];

MonotoneSolverFromData::usage =
"MonotoneSolverFromData[Data] solves the MFG problem from raw Data using the monotone operator method.
Options: \"TimeSteps\" (default 100) - number of ODE integration steps.";
MonotoneSolverFromData[Data_, opts:OptionsPattern[]] :=
	Module[{d2e},
		d2e = Data2Equations[Data];
		MonotoneSolver[d2e, opts]
	];

MonotoneSolver[d2e_, opts:OptionsPattern[]] :=
	Module[{B, K, cost, jj, cc, x0},
		{B, K, cost, jj} = GetKirchhoffMatrix[d2e];
		cc[j_?NumberVectorQ] := cost[j];
		x0 = jj /. (First@FindInstance[K . jj == B && And @@ ((# > 0) & /@ jj), jj]);
		MonotoneSolverODE[x0, K, jj, cc, GradientProjection, FilterRules[{opts}, Options[MonotoneSolverODE]]]
	];

Options[MonotoneSolverODE] = {"TimeSteps" -> 100};

MonotoneSolverODE::usage =
"MonotoneSolverODE[x0, K, jj, cc, gradProj] solves the gradient flow ODE from initial condition x0.
Options: \"TimeSteps\" (default 100).";
MonotoneSolverODE[x0_, K_, jj_, cc_, gradProj_, opts:OptionsPattern[]] :=
	Module[{At, dim, sol, x, xx, n},
		n = OptionValue["TimeSteps"];
		At = Transpose[K];
		dim = Length[x0];
		sol = NDSolve[
			x'[t] == -gradProj[x[t], K, dim, At] . cc[x[t]] && x[0] == x0,
			x[t], {t, 0, n}, Method -> "BDF"];
		xx = Table[x[t] /. sol[[1]], {t, 1, n, 1}];
		AssociationThread[jj, Last[xx]]
	];
