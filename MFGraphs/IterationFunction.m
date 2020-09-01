(* ::Package:: *)

(* Wolfram Language package *)


(*ToRules @ Or @@ Reduce[# && Data[type2] && Data[type3]]&/@Data[type1] (*if Data[type1] has head Or*)
*)

(*TODO learn how to document code!*)

(*TODO make the iterative function work on all the values.*)


(*TODO define usage properly*)
FixedSolverStepX2::usage = 
"";

EqEliminatorX::usage = 
"";

EqEliminatorX2::usage = 
"Improving the Eliminator.";

Begin["`Private`"]

EqEliminatorX[{system_, rules_}] := (*if it solves, system becomes True. rules is a list of rules or an association.*)
    Which[Head[system] === And, (*TODO what if the Head is Or? As in a boolean-converted system...*)
        Module[ {newrules, rulesAss = Association[rules]},(*replace rules in system before solving*)
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
        Print["EqEliminatorX: Reached a solution!"];
        {system, rules},
        True,
        Print["Eliminator: System is not True and Head of system is NOT And or Or. Returning the input."];
        {system, rules}
    ];
   (*TODO critical congestion solution tester: replace "solution" on the equations. Look at the output to figure out if it works!! *)
   (*in the nonlinear probelm, look for approximately equalities.*)
   (*TODO use the result of the fixed point in this function and the fixed point simplifier again.*)
  
EqEliminatorX2[{AA_, rules_}] :=
    Module[ {NN, OO, EE, AA2, Asrules = Association[rules]},
    	AA2 = AA/.Asrules;
        NN = Select[AA2, Head[#] === NonNegative &];
        OO = Select[AA2, Head[#] === Or &];
        EE = Select[AA2, Head[#] === Equal &];
        {OO, NN, EE, AssociateTo[Asrules, Solve[EE, Reals] // First // Quiet]}
    ]
  
EqEliminatorX2[{OO_, NN_, EE_, sol_Association}] :=
    Module[ {OO2, NN2, EE2, AA2, Assol = sol},
        NN2 = Reduce[(NN && OO)/. sol, Reals];
        (*Print["NonNegatives: ", NN2];*)
        OO2 = OO /. sol;
        AA2 = NN2 && OO2;
        If[AA2 === True, {OO2, NN2, EE /. Assol, AssociateTo[Assol, Solve[EE, Reals] // First // Quiet]/. Assol},
        	EE2 = Select[AA2, Head[#] === Equal &];
        	{OO2, NN2, EE2, AssociateTo[Assol, Solve[EE2, Reals] // First // Quiet]}
        ]
    ]

FixedSolverStepX2[Eqs_Association][rules_] :=
    Module[ {system, nonlinear, newsolve},
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], (Eqs["Nrhs"] /. rules)}]);
        (*Print["Step: nonlinear: ",nonlinear];(*hopefully these are all numeric.*)*)
        system = Eqs["reduced"][[1]] && nonlinear;
        (*Print["Step: SYSTEM: ", system];*)
        newsolve = FixedPoint[EqEliminatorX, {system, Eqs["reduced"][[2]]}, 10];
        (*Print["Step: newsolve 2: ", newsolve];*)
        (*Print["Step: Reducing ..."];(*be carefull here!*)(*TODO think of the cases when Reduce returns False...*)*)
        newsolve = {Reduce[newsolve[[1]],Reals],newsolve[[2]]}//Quiet;(*We use Quiet because of: Reduce::ratnz: Reduce was unable to solve the system with inexact coefficients. The answer was obtained by solving a corresponding exact system and numericizing the result.*)
        (*Print["Step: newsolve 2 resolved", newsolve];*)
        newsolve = FixedPoint[EqEliminatorX, newsolve, 10];
        (*Print["Step: newsolve 3: ", newsolve];*)
        Association @ newsolve[[2]]
    ]



End[]
