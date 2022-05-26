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
    EqPosJs, EqPosJts, Nrhs, Nlhs, style, styler, styleg},
   	 style = Style[#,Bold,Blue]&; 
     styler = Style[#,Bold,Red]&; 
     styleg = Style[#,Bold,Darker@Green]&; 
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
     bool = (EqEntryIn&&EqValueAuxiliaryEdges&&EqCompCon&&
     	EqBalanceSplittingCurrents&&EqCurrentCompCon&&EqTransitionCompCon&&
     	EqPosJs&&EqPosJts&&EqSwitchingByVertex)/.assoc;
     If[bool === True, Print["All restrictions are ", styleg@"True"],
     Print["At least one of the restrictions is ", styler@"False"];
     Print[style@"EqEntryIn: ", EqEntryIn/.assoc, "\n", EqEntryIn];
     Print[style@"EqValueAuxiliaryEdges: ", EqValueAuxiliaryEdges/.assoc, "\n", EqValueAuxiliaryEdges];
     Print[style@"EqCompCon: ", EqCompCon/.assoc, "\n", EqCompCon];
     Print[style@"EqBalanceSplittingCurrents: ", EqBalanceSplittingCurrents/.assoc, "\n", EqBalanceSplittingCurrents];
     Print[style@"EqCurrentCompCon: ", EqCurrentCompCon/.assoc, "\n",EqCurrentCompCon];
     Print[style@"EqTransitionCompCon: ", EqTransitionCompCon/.assoc, "\n", EqTransitionCompCon];
     Print[style@"EqPosJs: ", EqPosJs/.assoc, "\n", EqPosJs];
     Print[style@"EqPosJts: ", EqPosJts/.assoc, "\n", EqPosJts];
     Print[style@"EqSwitchingByVertex: ", EqSwitchingByVertex/.assoc, "\n", EqSwitchingByVertex]];
     Print[style@"Nrhs: ", N[Nrhs/.assoc], "\n", Nrhs];
     Print[style@"Nlhs: ", N[Nlhs/.assoc], "\n", Nlhs];
     Print[style@"Max error for non-linear solution: ", Norm[N[Nlhs/.assoc]-(Nrhs/.assoc),Infinity]];
     N@assoc
	];
	
RoundValues[x_?NumberQ] := Round[x, 10^-13]

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