(* Wolfram Language package *)

Hess::usage = 
"Hess[j] returns a numeric diagonal matrix with the reciprocals of the elements in the (numeric) vector j.";
Hess[j_?NumberVectorQ] := DiagonalMatrix[1/# & /@ j];

(*These, H2 and H3 are alternatives to the function Hess*)
H2[j_?NumberVectorQ] := DiagonalMatrix[1/#^2 & /@ j];

H3[j_?NumberVectorQ] := DiagonalMatrix[1/#^(1/4) & /@ j];


InvH::usage = 
"InvH[j] returns a diagonal matrix with the (numeric) vector j. This is the inverse of the matrix Hess[j].";
InvH[j_?NumberVectorQ] := DiagonalMatrix[j];

NumberMatrixQ::usage = 
"NumberMatrixQ[A] returns True if the elements of the matrix A are numeric.";
NumberMatrixQ[A_] := NumberVectorQ@Flatten@A;

InvAHA::usage = 
"InvAHA[j,A,At] returns the product of A times the inverse of Hess[j] times At";
InvAHA[x_?NumberVectorQ, A_, At_] := A  . InvH[x] . At;
InvAHA[x_?NumberVectorQ, A_] := A  . InvH[x] . Transpose[A];

GG::usage = (*this function can be optimized, If At is A transposed, and dim can be obtained from A*)
  "GG[x,A,dim,At] is a good matrix.";
GG[x_?NumberVectorQ, A_, dim_, At_] := 
  InvH[x] . (IdentityMatrix[dim] - At . PseudoInverse[InvAHA[x, A, At] ] . A . InvH[x] );(*Taking Pseudo, instead of just Inverse*)

GG[x_?NumberVectorQ, A_] :=
	Module[{dim, At},
		dim = Length[A];
		At = Transpose[A];
		InvH[x] . (IdentityMatrix[dim] - At . PseudoInverse[InvAHA[x, A, At] ] . A . InvH[x] );
	];
	
MonotoneSolverFromData::usage = 
"MonotoneSolverFromData[Data] ";
Options[MonotoneSolverData]={"n"->100};	
MonotoneSolverFromData[Data_]:= 
	Module[{d2e}, 
		d2e = Data2Equations[Data];
		MonotoneSolver[d2e]
	];
	
Options[MonotoneSolver]=Options[MonotoneSolverData];
MonotoneSolver[d2e_]:=
	Module[{B,K,cost,jj,cc,x0},
		{B,K,cost,jj} = GetKirchhoffMatrix[d2e];
		cc[j_?NumberVectorQ] := cost[j]; 
		x0 = jj /. (First@FindInstance[K . jj == B && And @@ ((# > 0) & /@ jj), jj]); (*Select a feasible initial condition for the current distribution*)
		MonotoneSolver[x0,K,jj,cc,GG]	
	];
MonotoneSolver[x0_,K_,jj_,cc_,GG_] := 
	Module[{At, dim, sol, x, xx, n},
		n=OptionValue["n"];
		At = Transpose[K];
		dim = Length[x0];
		sol = NDSolve[x'[t] == -GG[x[t], K, dim, At, x[t]].cc[x[t]] && x[0] == x0, x[t], {t, 0, n}, Method -> "BDF"];
		xx = Table[x[t] /. sol[[1]], {t, 1, n, 1}];
		AssociationThread[jj, xx[[n]]]
	];
	