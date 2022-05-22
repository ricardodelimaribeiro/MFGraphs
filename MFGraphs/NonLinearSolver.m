(*Wolfram Language package*)
Get["/Users/ribeirrd/eclipse-workspace/MFGraphs/MFGraphs/D2E2.m"];
NonLinear[Eqs_] := 
  Module[{PreEqs, Nlhs, Nrhs, EqApprox, InitRules, NewSystem, 
    approxrules}, PreEqs = MFGPreprocessing[Eqs];
   Print["Preprocessing done!"];
   Nlhs = Lookup[PreEqs, "Nlhs", $Failed];
   Nrhs = Lookup[PreEqs, "Nrhs", $Failed];
   approxrules = CriticalCongestionSolver[PreEqs];
   Print["raw: ", Nrhs];
   InitRules = Lookup[PreEqs, "InitRules", $Failed];
   Print["with initrules: ", Nrhs /. InitRules];
   Print["with j: ", Nrhs /. approxrules];
   EqApprox = And @@ MapThread[Equal, {Nlhs, Nrhs}];
   Print[EqApprox];
   EqApprox = And @@ MapThread[Equal, {Nlhs, Nrhs /. approxrules}];
   Print[EqApprox];
   NewSystem = Lookup[PreEqs, "NewSystem", $Failed];
   (*Updated (with rules from preprocessing) edge equations*)(*Print[
   "22: ",
   EqCritical];*){NewSystem, InitRules} = 
    TripleClean[{Join[{EqApprox}, Take[NewSystem, {2, 3}]], 
      InitRules}];
   (*Print["33: ",NewSystem];*)NewSystem = NewReduce[And @@ NewSystem];
   NewSystem = BooleanConvert[NewSystem, "CNF"];
   (*Print["44: ",NewSystem];*)
   NewSystem = 
    Sys2Triple[BooleanConvert[Reduce[NewSystem, Reals], "CNF"]];
   {NewSystem, InitRules} = TripleClean[{NewSystem, InitRules}];
   If[And @@ NewSystem === False, Print["There is no solution"]];
   InitRules];

NonLinearStep[Eqs_][approxrules_] := 
  Module[{Nlhs, EqNonLinear, InitRules, NewSystem}, 
   Nlhs = Lookup[Eqs, "Nlhs", $Failed];
   InitRules = Lookup[Eqs, "InitRules", $Failed];
   NewSystem = Lookup[Eqs, "NewSystem", $Failed];
   InitRules/.approxrules...
   (*Updated (with rules from preprocessing) edge equations*)
   EqNonLinear = (And @@ ((# == 0) & /@ Nlhs)) /. InitRules;
   (*Print["FinalClean..."];*)
   MFGSystemSolver[Eqs][EqNonLinear]
   ];
   
   

RoundValues[x_?NumberQ] := Round[x, 10^-10]

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