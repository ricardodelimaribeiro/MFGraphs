(*Wolfram Language package*)
Get["/Users/ribeirrd/eclipse-workspace/MFGraphs/MFGraphs/D2E2.m"];
NonLinear[Eqs_] :=    
  Module[{AssoCritical, PreEqs, AssoNonCritical, js, InitRules, MaxIter = 50}, 
   PreEqs = If[KeyExistsQ[Eqs,"InitRules"], Eqs, MFGPreprocessing[Eqs]];
   AssoCritical = Lookup[Eqs, "AssoCritical", CriticalCongestionSolver[Eqs]];
   js = Lookup[Eqs, "js",$Failed];
   AssoNonCritical = FixedPointList[NonLinearStep[PreEqs], KeyTake[AssoCritical,js], MaxIter];
   Print["Iterated ", Length[AssoNonCritical], " times out of ", MaxIter];
   AssoNonCritical = AssoNonCritical//Last;
   InitRules = Lookup[Eqs, "InitRules", $Failed];
   AssoNonCritical = MFGSystemSolver[Eqs][AssoNonCritical];
   (*Print[AssoNonCritical];*)
   (*costpluscurrents = Lookup[Eqs, "costpluscurrents",$Failed];
   Print[InitRules];
   Print[InitRules/.costpluscurrents/.AssoNonCritical];*)
   Join[PreEqs, Association["AssoNonCritical" -> AssoNonCritical]]
   ];

NonLinearStep[Eqs_][approxJs_] := 
	Module[{approx, js},
		(*Print[approxJs];*)
		js = Lookup[Eqs, "js",$Failed];
		approx = MFGSystemSolver[Eqs][approxJs];
		approx = KeyTake[approx, js];
		(*Print["step: ",approx];*)
		approx
	];
	(*TODO: assoc is given the right way? can we get it from the Eqs?*)
IsNonLinearSolution[Eqs_][assoc_]:=
Module[{EqEntryIn, EqValueAuxiliaryEdges, EqSwitchingByVertex, EqCompCon, 
    EqBalanceSplittingCurrents, EqCurrentCompCon, EqTransitionCompCon,
    EqPosJs, EqPosJts}, 
     EqEntryIn = And @@ Lookup[Eqs, "EqEntryIn", $Failed];
     EqCompCon = Lookup[Eqs, "EqCompCon", $Failed];
     EqCurrentCompCon = Lookup[Eqs, "EqCurrentCompCon", $Failed];
     EqTransitionCompCon = Lookup[Eqs, "EqTransitionCompCon", $Failed];
     EqPosJs = Lookup[Eqs, "EqPosJs", $Failed];
     EqPosJts = Lookup[Eqs, "EqPosJts", $Failed];
     EqSwitchingByVertex = And @@ Lookup[Eqs, "EqSwitchingByVertex", $Failed];
     EqBalanceSplittingCurrents = Lookup[Eqs, "EqBalanceSplittingCurrents", $Failed];
     EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", $Failed];
     EqGeneral = Lookup[Eqs, "EqGeneral", $Failed];
     Nrhs = Lookup[Eqs, "Nrhs", $Failed];
     Nlhs = Lookup[Eqs, "Nlhs", $Failed];
     costpluscurrents = Lookup[Eqs, "costpluscurrents", $Failed];
     js = Lookup[Eqs, "js", $Failed];
     
     Print["EqEntryIn: ", EqEntryIn, "\t", EqEntryIn/.assoc];
     Print["EqValueAuxiliaryEdges: ", EqValueAuxiliaryEdges, "\t",EqValueAuxiliaryEdges/.assoc];
     Print["EqCompCon: ", EqCompCon, "\t", EqCompCon/.assoc];
     Print["EqBalanceSplittingCurrents: ", EqBalanceSplittingCurrents, "\t",EqBalanceSplittingCurrents/.assoc];
     Print["EqCurrentCompCon: ", EqCurrentCompCon, "\t",EqCurrentCompCon/.assoc];
     Print["EqTransitionCompCon: ", EqTransitionCompCon, "\t",EqTransitionCompCon/.assoc];
     Print["EqPosJs: ", EqPosJs, "\t",EqPosJs/.assoc];
     Print["EqPosJts: ", EqPosJts, "\t",EqPosJts/.assoc];
     Print["EqSwitchingByVertex: ", EqSwitchingByVertex, "\t",Simplify/@(EqSwitchingByVertex/.assoc)];
     Print["Nlhs: ", Nlhs, "\t", N[Nlhs/.assoc]];
     Print["Nrhs: ", Nrhs, "\t", Nrhs/.assoc];
     Print["Max error for non-linear solution: ", Max[N[Nlhs/.assoc]-(Nrhs/.assoc)]];
     assoc
	];
	
RoundValues[x_?NumberQ] := Round[x, 10^-11]

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