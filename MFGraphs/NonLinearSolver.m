(* Wolfram Language package *)
Get["/Users/ribeirrd/eclipse-workspace/MFGraphs/MFGraphs/D2E2.m"];
NonLinear[Eqs_] := 
  Module[{PreEqs, Nlhs, Nrhs, EqCritical, InitRules, NewSystem}, 
   PreEqs = MFGPreprocessing[Eqs];
   Print["Preprocessing done!"];
   Nlhs = Lookup[PreEqs, "Nlhs", $Failed];
   Nrhs = Lookup[PreEqs, "Nrhs", $Failed];
   Print[Nrhs];
   InitRules = Lookup[PreEqs, "InitRules", $Failed];
   Print[Nrhs/.InitRules];
   NewSystem = Lookup[PreEqs, "NewSystem", $Failed];
   (*Updated (with rules from preprocessing) edge equations*)
   EqCritical = (And @@ ((# == 0) & /@ Nlhs)) /. InitRules;
   (*Print["22: ",
   EqCritical];*){NewSystem, InitRules} = 
    TripleClean[{Join[{EqCritical}, Take[NewSystem, {2, 3}]], 
      InitRules}];
   (*Print["33: ",NewSystem];*)NewSystem = NewReduce[And @@ NewSystem];
   NewSystem = BooleanConvert[NewSystem, "CNF"];
   (*Print["44: ",NewSystem];*)
   NewSystem = 
    Sys2Triple[BooleanConvert[Reduce[NewSystem, Reals], "CNF"]];
   {NewSystem, InitRules} = TripleClean[{NewSystem, InitRules}];
   If[And @@ NewSystem === False, Print["There is no solution"]];
   InitRules];

RoundValues[x_?NumberQ] :=
    Round[x, 10^-10]
    
RoundValues[Rule[a_, b_]] :=
    Rule[a, RoundValues[b]]
    
RoundValues[x_List] :=
    RoundValues /@ x
    
RoundValues[x_Association] :=
    RoundValues /@ x
