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

EqEliminatorX[system_] :=
Module[ {newrules},(*replace rules in system before solving*)
 	newrules = Select[system, (Head[#] === Equal) &] && And @@ (globalrules /. Rule -> Equal) // Solve // First // Quiet;
    (*Print[newrules]*);
    globalrules =  newrules;(*replace in the rules before appending*)
    system /. globalrules
];
   
StartSolverX[Eqs_Association] :=
Module[ {eqs, listeqs, rules},
	eqs = DeleteDuplicates[Select[Reduce[# && And @@ ((# == 0) & /@ Eqs["Nlhs"] /. Eqs["globalrules"]), Reals, Backsubstitution -> True] & /@ Eqs["EqAllAll"], (# =!= False) &]];
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
	nonlinear = And @@ (MapThread[(#1 == #2) &, {Eqs["Nlhs"] /. Eqs["globalrules"], (Eqs["Nrhs"] /. rules[[1]])}]);
    Print[nonlinear];
    system = nonlinear && Eqs["reduced"];
    Print[system];
    globalrules = Eqs["globalrules"];
    newsolve = FixedPoint[EqEliminatorX, First@system, 10];
    Print[newsolve];
    globalrules
       (*Try to use EqEliminatorX here. 
           But rembemer it changes the globalrules variable (which is global)*)
]
  
End[]
