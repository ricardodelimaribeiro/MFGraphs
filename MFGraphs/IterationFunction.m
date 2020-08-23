(* ::Package:: *)

(* Wolfram Language package *)


(*ToRules @ Or @@ Reduce[# && Data[type2] && Data[type3]]&/@Data[type1] (*if Data[type1] has head Or*)
*)

(*TODO learn how to document code!*)

(*TODO make the iterative function work on all the values.*)

EqEliminatorX::usage = 
"";
StartSolverX::usage = 
"";
FixedPointSolverStepX::usage = 
"";

Begin["`Private`"]

EqEliminatorX[{system_, rules_}] :=
If[Head[system] === And, 
	Module[ {newrules},(*replace rules in system before solving*)
 	newrules = Select[system, (Head[#] === Equal) &] && And @@ (rules /. Rule -> Equal) // Solve // First // Quiet;
    {system /. newrules, newrules}],
    {system, rules}
];
   (*TODO critical congestion solution tester: replace "solution" on the equations. Look at the output to figure out if it works!! *)
   (*in the nonlinear probelm, look for approximately equalities.*)
   (*TODO use the result of the fixed point in this function and the fixed point simplifier again.*)
StartSolverX[Eqs_Association] :=
Module[ {eqs, listeqs, rules},
	eqs = DeleteDuplicates[Select[Reduce[# && And @@ ((# == 0) & /@ Eqs["Nlhs"] /. Eqs["reducing rules"]), Reals, Backsubstitution -> True] & /@ Eqs["EqAllAll"], (# =!= False) &]];
    listeqs = List @@ eqs;
    rules = Solve[List @@ eqs, Reals];
    If[ Length[rules] > 1,
    	Print["There is more than one set of NotFalse equations!"];
        rules,
        rules
    ]
]
  
FixedPointSolverStepX[Eqs_Association][rules_] :=
Module[ {system, nonlinear, newsolve},
	nonlinear = And @@ (MapThread[(#1 == #2) &, {Eqs["Nlhs"] /. Eqs["reducing rules"], (Eqs["Nrhs"] /. rules[[1]])}]);
    Print[nonlinear];
    system = nonlinear && Eqs["reduced system"];
    Print[system];
    globalrules = Eqs["reducing rules"];
    newsolve = FixedPoint[EqEliminatorX, {system, rules}, 10];
    
    	
    Print[newsolve];
    globalrules
       (*Try to use EqEliminatorX here. 
           But rembemer it changes the globalrules variable (which is global)*)
]
  
End[]
