(* ::Package:: *)

(* Wolfram Language package *)


(*ToRules @ Or @@ Reduce[# && Data[type2] && Data[type3]]&/@Data[type1] (*if Data[type1] has head Or*)
*)

(*TODO learn how to document code!*)

(*TODO make the iterative function work on all the values.*)

FixedSolverStepX2::usage = 
"";

EqEliminatorX::usage = 
"";

StartSolverX::usage = 
"";

FixedPointSolverStepX::usage = 
"";

Begin["`Private`"]

EqEliminatorX[{system_, rules_}] := (*if it solves, system becomes True. rules is a list of rules or an association.*)
Which[Head[system] === And, 
	Module[ {newrules, rulesAss},(*replace rules in system before solving*)
		rulesAss = Association[rules];
		newrules = Select[system, (Head[#] === Equal) &]  // Solve // First // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)
		(*Print["newrules :", newrules]*);
    	{system /. newrules, AssociateTo[rulesAss, newrules]/.newrules}],
    system === True,
    	Print["EqEliminatorX: Reached a solution!"];
    	{system, rules},
    True,
    	Print["Something is wrong with the equations: ", system];    
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
    	Print["There is more than one set of NotFalse equations!"];(*rules is a list with solution lists*)
        rules,
        rules
    ]
]
  
FixedPointSolverStepX[Eqs_Association][rules_] :=
Module[ {system, nonlinear, newsolve, globalrules},
	nonlinear = And @@ (MapThread[(#1 == #2) &, {Eqs["Nlhs"] /. Eqs["reducing rules"], (Eqs["Nrhs"] /. rules)}]);
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
  
FixedSolverStepX2[Eqs_Association][rules_] :=
Module[ {system, nonlinear, newsolve},
	nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], (Eqs["Nrhs"] /. rules)}]);
    (*Print["Step: nonlinear: ",nonlinear];(*hopefully these are all numeric.*)*)
    system = Eqs["reduced"][[1]] && nonlinear;
    (*Print["Step: SYSTEM: ", system];*)
    newsolve = FixedPoint[EqEliminatorX, {system, Eqs["reduced"][[2]]}, 10];
    (*Print["Step: newsolve 2: ", newsolve];*)
    Print["Step: Reducing ..."];
    newsolve = {newsolve[[1]]//Reduce,newsolve[[2]]}//Quiet;(*We use Quiet because of: Reduce::ratnz: Reduce was unable to solve the system with inexact coefficients. The answer was obtained by solving a corresponding exact system and numericizing the result.*)
    newsolve = FixedPoint[EqEliminatorX, newsolve,10];
    (*Print["Step: newsolve 3: ", newsolve];*)
	Association @ newsolve[[2]]
]



End[]
