(*Wolfram Language package*)
Get["/Users/ribeirrd/eclipse-workspace/MFGraphs/MFGraphs/D2E2.m"];
(*NonLinear[Eqs_][initial_] :=
	Module[{AssoNonCritical, MaxIter = 50},
		AssoNonCritical = FixedPointList[NonLinearStep[Eqs], initial, MaxIter];
		Print["Iterated ", Length[AssoNonCritical], " times out of ", MaxIter];
		AssoNonCritical
	];
*)
NonLinear::usage =
	"NonLinear[Eqs] takes an association ";

NonLinear[Eqs_] :=    
  Module[{AssoCritical, PreEqs = Eqs, AssoNonCritical, NonCriticalList, js = 1, MaxIter = 15}, 
   If[KeyExistsQ[PreEqs, "AssoCritical"], 
   	(*if there is already an approximation for the non congestion case, use it!*)
   	AssoCritical = Lookup[PreEqs, "AssoNonCritical", PreEqs["AssoCritical"]],
   	PreEqs = MFGPreprocessing[PreEqs];
   	js = Lookup[PreEqs, "js", $Failed];
   	AssoCritical = AssociationThread[js, 0 js]
   ];
   NonCriticalList = FixedPointList[NonLinearStep[PreEqs], AssoCritical, MaxIter];
   Print["Iterated ", Length[NonCriticalList]-1, " times out of ", MaxIter];
   If[js =!= 1, AssoCritical = NonCriticalList[[2]]];
   AssoNonCritical = NonCriticalList // Last;
   Join[PreEqs, Association[{"AssoCritical" -> AssoCritical, "AssoNonCritical" -> AssoNonCritical}]]
   ];

NonLinearStep[Eqs_][approxSol_] := 
	Module[{approxJs, approx, js, Nrhs, Nlhs},
		js = Lookup[Eqs, "js", $Failed];
		Nrhs = Lookup[Eqs, "Nrhs", $Failed];
     	Nlhs = Lookup[Eqs, "Nlhs", $Failed];
     	approxJs = KeyTake[approxSol, js];
		approx = MFGSystemSolver[Eqs][approxJs];
		Print["Max error for non-linear solution: ", Norm[N[Nlhs/.approx]-(Nrhs/.approx),Infinity]];
		approx
	];
	
	(*TODO: assoc is given the right way? can we get it from the Eqs?*)
IsNonLinearSolution[Eqs_][assoc_]:=
Module[{EqEntryIn, EqValueAuxiliaryEdges, EqSwitchingByVertex, EqCompCon, 
    EqBalanceSplittingCurrents, EqCurrentCompCon, EqTransitionCompCon,
    EqPosJs, EqPosJts, Nrhs, Nlhs}, 
     EqEntryIn = And @@ Lookup[Eqs, "EqEntryIn", $Failed];
     EqCompCon = Lookup[Eqs, "EqCompCon", $Failed];
     EqCurrentCompCon = Lookup[Eqs, "EqCurrentCompCon", $Failed];
     EqTransitionCompCon = Lookup[Eqs, "EqTransitionCompCon", $Failed];
     EqPosJs = Lookup[Eqs, "EqPosJs", $Failed];
     EqPosJts = Lookup[Eqs, "EqPosJts", $Failed];
     EqSwitchingByVertex = And @@ Lookup[Eqs, "EqSwitchingByVertex", $Failed];
     EqBalanceSplittingCurrents = Lookup[Eqs, "EqBalanceSplittingCurrents", $Failed];
     EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", $Failed];
     Nrhs = Lookup[Eqs, "Nrhs", $Failed];
     Nlhs = Lookup[Eqs, "Nlhs", $Failed];
     Print["All restrictions: ", (EqEntryIn&&EqValueAuxiliaryEdges&&EqCompCon&&
     	EqBalanceSplittingCurrents&&EqCurrentCompCon&&EqTransitionCompCon&&
     	EqPosJs&&EqPosJts&&EqSwitchingByVertex)/.assoc];
     Print["EqEntryIn: ", EqEntryIn/.assoc, "\n", EqEntryIn];
     Print["EqValueAuxiliaryEdges: ", EqValueAuxiliaryEdges/.assoc, "\n", EqValueAuxiliaryEdges];
     Print["EqCompCon: ", EqCompCon/.assoc, "\n", EqCompCon];
     Print["EqBalanceSplittingCurrents: ", EqBalanceSplittingCurrents/.assoc, "\n", EqBalanceSplittingCurrents];
     Print["EqCurrentCompCon: ", EqCurrentCompCon/.assoc, "\n",EqCurrentCompCon];
     Print["EqTransitionCompCon: ", EqTransitionCompCon/.assoc, "\n", EqTransitionCompCon];
     Print["EqPosJs: ", EqPosJs/.assoc, "\n", EqPosJs];
     Print["EqPosJts: ", EqPosJts/.assoc, "\n", EqPosJts];
     Print["EqSwitchingByVertex: ", EqSwitchingByVertex/.assoc, "\n", EqSwitchingByVertex];
     Print["Nlhs: ", N[Nlhs/.assoc], "\n", Nlhs];
     Print["Nrhs: ", Nrhs/.assoc, "\n", Nrhs];
     Print["Max error for non-linear solution: ", Norm[N[Nlhs/.assoc]-(Nrhs/.assoc),Infinity]];
     assoc
	];
	
RoundValues[x_?NumberQ] := Round[x, 10^-15]

RoundValues[Rule[a_, b_]] := Rule[a, RoundValues[b]]

RoundValues[x_List] := RoundValues /@ x

RoundValues[x_Association] := RoundValues /@ x

alpha[edge_]:= 1.2

g[m_, edge_]:= m

A = 0;

W::usage =
""
W = Function[{y, a}, a Sin[2 Pi (y + 1/4)]^2];

V::usage = 
""
V = Function[{x, edge}, W[x, A]];

H::usage = (*TODO include the other parameter here, if we call it beta, change the the other beta!*)
""
(*H = Function[{xi,p,m, edge}, p^2/(2 m^alpha[edge]) + V[xi, edge] - g[m]];*)
H = Function[{xi,p,m, edge}, p^2/(2 m^alpha[edge]) + V[xi, edge] - g[m, edge]];


M[j_?NumericQ, x_?NumericQ, edge_] := 
 If[j == 0. || j == 0, 0., 
  Values@First@FindRoot[H[x, -j/(m^(1 - alpha[edge])), m, edge], {m, 1}]]

IntM[j_?NumericQ, edge_] := 
 If[j == 0. || j == 0, 0., 
  j NIntegrate[M[j, x, edge]^(alpha[edge] - 1), {x, 0, 1}] // Quiet]