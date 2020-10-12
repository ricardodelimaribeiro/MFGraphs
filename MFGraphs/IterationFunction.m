(* ::Package:: *)

(* Wolfram Language package *)


(*ToRules @ Or @@ Reduce[# && Data[type2] && Data[type3]]&/@Data[type1] (*if Data[type1] has head Or*)
*)

(*TODO learn how to document code!*)

(*TODO make the iterative function work on all the values.*)


(*TODO define usage properly*)
FixedSolverStepX2::usage = 
"";

FixedSolverStepX1::usage = 
"";

EqEliminatorX::usage = 
"";

EqEliminatorX2::usage = 
"Improving the Eliminator.";

Begin["`Private`"]

EqEliminatorX[{system_, rules_}] := (*if it solves, system becomes True. rules is a list of rules or an association.*)
	(*system = system /. rules;*)(*TODO why does this not work?*)
	(*Print[system];*)
    Which[Head[system] === And, (*TODO what if the Head is Or? As in a boolean-converted system...*)
        Module[ {newrules, rulesAss = Association[rules]},(*replace rules in system before solving*)
        	(*system = system /. rules;*)(*TODO why not use rules in the system before anything else?*)
            newrules = Select[system, (Head[#] === Equal) &]  // Solve // First // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)
            (*Print["rulesAss: ", rulesAss];
            Print["newrules: ", newrules];
            Print["system: ",system];
            Print["system newrules: ", system /. newrules];*)
            {system /. newrules, AssociateTo[rulesAss, newrules]/.newrules}
        ],
        Head[system] === Equal,
        Module[ {newrules, rulesAss = Association[rules]},(*replace rules in system before solving*)
            newrules = First[Solve[system]];
            {system /. newrules, AssociateTo[rulesAss, newrules]/.newrules}
        ],        
        system === True,
        (*Print["EqEliminatorX: Reached a solution!"];
        *){system, rules},
        True,
        Print["Eliminator: System is not True and Head of system is NOT And or Equal. Returning the input."];
        {system, rules}
    ];
   (*TODO critical congestion solution tester: replace "solution" on the equations. Look at the output to figure out if it works!! *)
   (*in the nonlinear probelm, look for approximately equalities.*)
   (*TODO use the result of the fixed point in this function and the fixed point simplifier again.*)
  
EqEliminatorX2[{AA_, rules_}] :=
    Module[ {NN, OO, EE, AA2 = AA, Asrules = Association[rules]},
    	(*Print["System Heads: ",List @@ DeleteDuplicates@(Head[#] & /@ AA)];
    	Print["Rules association input: ", Asrules];*)
    	(*AA2 = AA /. Asrules;*)
    	EE = Select[AA2, Head[#] === Equal &];
        NN = Select[AA2, Head[#] === NonNegative &];
        OO = Select[AA2, Head[#] === Or &];
        (*Print["NonNegative: ", NN];
        Print["Alternatives: ", OO];
        Print["Equalities: ", EE];*)
        AssociateTo[Asrules, Solve[EE, Reals] // First // Quiet];
        (*Print[Asrules];*)
        {OO, NN, Asrules/. Asrules} 
    ]

EqEliminatorX2[{OO_, NN_, sol_Association}] :=
    Module[ {OO2, NN2, EE2, AA2, Assol = sol},
        NN2 = NN /. Assol;
        OO2 = OO /. Assol;
        AA2 = BooleanConvert[Reduce[NN2 && OO2, Reals] // Quiet,"CNF"]; (*Quiet : Reduce::ratnz: Reduce was unable to solve the system with inexact coefficients. The answer was obtained by solving a corresponding exact system and numericizing the result.*)
        (*Print["System Heads: ",List @@ DeleteDuplicates@(Head[#] & /@ AA2)];*)
    	(*Print["Rules association input: ", Assol];*)
        (*Print["System Heads: ",List @@ DeleteDuplicates@(Head[#] & /@ AA2)];
        Print["NonNegative: ", NN2];
        Print["Alternatives: ", OO2];*)
        (*Print["EqEliminatorX2: Sistem: ", AA2];*)
        If[AA2 === True, 
        	{OO2, NN2, Assol},
        	EE2 = Select[AA2, Head[#] === Equal &];
        	AssociateTo[Assol, Solve[EE2, Reals] // First // Quiet];
        	{OO2, NN2, Assol} /. Assol
        ]
    ]
    
 

FixedSolverStepX2[Eqs_Association][rules_] :=
    Module[ {system, nonlinear, newsolve},
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[(Eqs["Nrhs"] /. rules),10^-14]}]);
        (*Print["FixedSolveX2: Step: nonlinear: ",nonlinear];(*hopefully these are all numeric.*)*)
        system = Eqs["EqAllAll"] && nonlinear;
        (*Print["Step: SYSTEM: ", system];*)
        newsolve = FixedPoint[EqEliminatorX2, {system, Eqs["reduced"][[2]]}, 10];
        (*Print["Step: newsolve 2: ", newsolve];*)
        (*Print["Step: Reducing ..."];(*be carefull here!*)(*TODO think of the cases when Reduce returns False...*)*)
        (*Print["Step: newsolve 3: ", newsolve];*)
        newsolve[[3]]
    ]
    
    (*TODO use Chop with a predefined tolerance*)
(*This is the implementation we were using that works in some cases.*)
FixedSolverStepX1[Eqs_Association][rules_] :=
    Module[ {system, nonlinear, newsolve},
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[(Eqs["Nrhs"] /. rules),10^-14]}]);
        (*Print["Step: nonlinear: ",nonlinear];(*hopefully these are all numeric.*)*)
        system = Eqs["EqAllAll"] && nonlinear;
        (*Print["Step: SYSTEM: ", system];*)
        newsolve = FixedPoint[EqEliminatorX, {system, Eqs["reduced"][[2]]}, 10];
        (*Print["Step: newsolve 2: ", newsolve];*)
        (*Print["Step: Reducing ..."];(*be carefull here!*)(*TODO think of the cases when Reduce returns False...*)*)
        newsolve = {Reduce[newsolve[[1]],Reals],newsolve[[2]]}//Quiet;(*We use Quiet because of: Reduce::ratnz: Reduce was unable to solve the system with inexact coefficients. The answer was obtained by solving a corresponding exact system and numericizing the result.*)
        (*Print["Step: newsolve 2 resolved ", newsolve];*)
        newsolve = FixedPoint[EqEliminatorX, newsolve, 10];
        (*Print["Step: newsolve 3: ", newsolve];*)
        Association @ newsolve[[2]]
    ]


End[]
