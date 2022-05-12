(* Wolfram Language package *)
Get["/Users/ribeirrd/eclipse-workspace/MFGraphs/MFGraphs/D2E2.m"];
NonLinear[Eqs_][approxrules_] := 
  Module[{PreEqs, Nlhs, Nrhs, EqApprox, InitRules, NewSystem}, 
   PreEqs = MFGPreprocessing[Eqs];
   Print["Preprocessing done!"];
   Nlhs = Lookup[PreEqs, "Nlhs", $Failed];	
   Nrhs = Lookup[PreEqs, "Nrhs", $Failed];
   Print["raw: ", Nrhs];
   InitRules = Lookup[PreEqs, "InitRules", $Failed];
   Print["with initrules: ",Nrhs /. InitRules];
   Print["with j: ", Nrhs/.approxrules];
   EqApprox = And @@ MapThread[Equal, {Nlhs,Nrhs}];
   Print[EqApprox];
   EqApprox = And @@ MapThread[Equal, {Nlhs,Nrhs/.approxrules }];
   Print[EqApprox];
   NewSystem = Lookup[PreEqs, "NewSystem", $Failed];
   (*Updated (with rules from preprocessing) edge equations*)
   (*Print["22: ",
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

NonLinearStep[Eqs_][j_]:=   
Module[{PreEqs, Nlhs, EqCritical, InitRules, NewSystem}, 
   PreEqs = MFGPreprocessing[Eqs];
   Nlhs = Lookup[PreEqs, "Nlhs", $Failed];
   InitRules = Lookup[PreEqs, "InitRules", $Failed];
   NewSystem = Lookup[PreEqs, "NewSystem", $Failed];
   (*Updated (with rules from preprocessing) edge equations*)
   EqCritical = (And @@ ((# == 0) & /@ Nlhs)) /. InitRules;
   Print["FinalClean..."];
   {NewSystem, InitRules} = 
   FinalClean[{Join[{EqCritical}, Take[NewSystem, {2, 3}]], 
     InitRules}];
   (*Print["33: ", {NewSystem, InitRules}];*)
   (*NewSystem = NewReduce[And @@ NewSystem];
   Print["34: ", NewSystem];
   NewSystem = BooleanConvert[NewSystem, "CNF"];
   
   NewSystem = 
    Sys2Triple[BooleanConvert[Reduce[NewSystem, Reals], "CNF"]];
   Print["55: ",NewSystem];
   {NewSystem, InitRules} = TripleClean[{NewSystem, InitRules}];
   Print["66: ",NewSystem];
   {NewSystem, InitRules} = TripleClean[{Sys2Triple[Reduce[NewSystem,Reals]],InitRules}];*)
   NewSystem = And @@ NewSystem;
   Which[NewSystem === False, Print["There is no solution"],
   	NewSystem =!= True, Print["Multiple solutions: ", NewSystem]
   ];
   InitRules];   

RoundValues[x_?NumberQ] :=
    Round[x, 10^-10]
    
RoundValues[Rule[a_, b_]] :=
    Rule[a, RoundValues[b]]
    
RoundValues[x_List] :=
    RoundValues /@ x
    
RoundValues[x_Association] :=
    RoundValues /@ x

M[j_?NumericQ, x_?NumericQ, edge_] :=
    If[ j == 0.|| j == 0,
        0.,
        Values @ First @ FindRoot[H[x, -j/(m^(1 - alpha)),m, edge], {m, 1}]
    ]
    
    IntM[j_?NumericQ, edge_] :=
    If[ j == 0.|| j == 0,
        0.,
        j NIntegrate[ M[j, x, edge]^(alpha-1), {x, 0, 1}] // Quiet
    ]