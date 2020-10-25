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

FixedSolverStepX3::usage = 
"";

FixedSolverStepX4::usage = 
"";

EqEliminatorX::usage = 
"";

EqEliminatorX2::usage = 
"Improving the Eliminator.";

EqEliminatorX3::usage = 
"Improving the Eliminator.";

Begin["`Private`"]

   (*TODO critical congestion solution tester: replace "solution" on the equations. Look at the output to figure out if it works!! *)
   (*in the nonlinear probelm, look for approximately equalities.*)
   (*TODO use the result of the fixed point in this function and the fixed point simplifier again.*)
  
EqEliminatorX2[{AA_, rules_List}] :=
    Module[ {EE, AA2, newrules, Asrules = Association @ rules, number},
    	(*Print["rules List"];*)
    	(*Print["System Heads: ",List @@ DeleteDuplicates@(Head[#] & /@ AA)];*)
    	EE = Select[AA, Head[#] === Equal &];
        AA2 = Select[AA, Head[#] =!= Equal &];
        (*Print["Equalities: ", EE];*)
        newrules = Solve[EE, Reals]//Quiet;
        number = Length @ newrules;
        If[number > 1,
        	Print["EqEliminatorX2: Number of solutions: ", Length @ newrules]];
        AssociateTo[Asrules, newrules // First];
        {AA2, Asrules}/. Asrules 
    ]
    
EqEliminatorX2[{AA_, rules_Association}] := 
Module[{EE, AA2, newrules, Asrules = rules, number},
		(*Print["AA: ", AA];*)
        AA2 = BooleanConvert[Reduce[AA, Reals] // Quiet,"CNF"]; (*Quiet : Reduce::ratnz: Reduce was unable to solve the system with inexact coefficients. The answer was obtained by solving a corresponding exact system and numericizing the result.*)
    	Print["EqEliminartorX2: Count! "];
        Which[AA2 === True, 
    		{AA2, Asrules},
        	AA2 === False,
        	Print["The system is False. Returning the inputs."];
        	{AA2, Asrules},
        	True,
        	EE = Select[AA2, Head[#] === Equal &];(*TODO what if there is no equality in AA2?*)
        	If[EE === True,
        		{AA2, Asrules},
        		newrules = Solve[EE, Reals]//Quiet;
        		(*Print[newrules];*)
        		number = Length[newrules];
        		Which[number > 1,
        			Print["EqEliminatorX2: Number of solutions: ", Length @ newrules, "\nThe system is: ", AA2],
        			number < 1,
        			Print["No solution for the system: ", EE];
        		];
        		AssociateTo[Asrules, newrules // First // Quiet];
        		(*Print[Asrules];*)
        		{AA2, Asrules} /. Asrules (*TODO return message and STOP if Total @ Values[Asrules] is not numeric*)
        	]
        ]
    ]
 
(*EqEliminatorX2[{AA_, rules_}] := 
    Module[ {EE},
        EE = Select[AA, Head[#] === Equal &];
        If[EE === True,
        	Print["There are no equalities."];
        	EqEliminatorX2[{AA,Association @ rules}],
        	Print["There are equalities!!"];
        	EqEliminatorX2[{AA,Normal @ rules}]
        ]
    ]*)
 

FixedSolverStepX2[Eqs_Association][rules_] :=
    Module[ {system, nonlinear, newsolve},
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[(Eqs["Nrhs"] /. rules),10^-6]}]);
        system = Eqs["EqAllAll"] && nonlinear;
        (*give system with sol to structural equations to equalities eliminator*)
        newsolve = FixedPoint[EqEliminatorX2, {system, Eqs["BoundaryRules"]}];
(*        Print["FixedX2: newsolve: ", {newsolve[[1]], Values@KeySort@newsolve[[2]]}];*)
        Print["FixedX2: next step: ", And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[(Eqs["Nrhs"] /. newsolve[[2]]),10^-6]}])];
        newsolve[[2]]
    ]
    

(*EqEliminatorX[{system_, rules_}] := (*if it solves, system becomes True. rules is a list of rules or an association.*)
Module[{newrules, rulesAss = Association @ rules},
	Which[Head[system] === And, (*TODO what if the Head is Or? As in a boolean-converted system...*)
		(*system = system /. rules;*)(*TODO why not use rules in the system before anything else?*)
			newrules = Select[system, (Head[#] === Equal) &]  // Solve // First // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)
        	(*Print["rulesAss: ", rulesAss];
        	Print["newrules: ", newrules];
        	Print["system: ",system];
        	Print["system newrules: ", system /. newrules];*)
        	{system, AssociateTo[rulesAss, newrules]} /. newrules,
        Head[system] === Equal,
        	newrules = First[Solve[system]];
        	{system, AssociateTo[rulesAss, newrules]} /. newrules,        
        system === True,
        (*Print["EqEliminatorX: Reached a solution!"];*)
        	{system, rules},
        True,
        	Print["Eliminator: System is not True and Head of system is NOT And or Equal. Returning the input."];
        	{system, rules}
    	]
	]*)

EqEliminatorX[{system_, rules_}] := (*if it solves, system becomes True. rules is a list of rules or an association.*)
    Which[
    	Head[system] === And, 
    		(*Print["And"];*)
        	Module[ {EE, newrules, rulesAss = Association[rules]},
        		EE = Select[system, (Head[#] === Equal) &];
            	If[EE === {}, 
            		{system, rules},
               		newrules =  Solve[EE, Reals] // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)
               		If[newrules === {},
               			{system, rules},
               			newrules = First @ newrules;  
               			{system, AssociateTo[rulesAss, newrules]}/.newrules
               		] 
            	]
        	],
        Head[system] === Equal,
        	Print["Equal"];
        	Module[ {newrules, rulesAss = Association[rules]},
            	newrules = Solve[system, Reals];
            	Print[newrules];
            	newrules = First @ newrules;
            	{system, AssociateTo[rulesAss, newrules]}/.newrules
        	],        
        system === True,
        	{system, rules},
        True,
        	Print["Eliminator: System is not True and Head of system is NOT And or Equal. Returning the input."];
        	{system, rules}
    ];
    
    
    (*TODO use Chop with a predefined tolerance*)
FixedSolverStepX1[Eqs_Association][rules_] :=
    Module[ {system, nonlinear, sol},
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[(Eqs["Nrhs"] /. rules),10^-6]}]);
        system = Eqs["EqAllAll"] && nonlinear;
        {system, sol} = FixedPoint[EqEliminatorX, {system, Eqs["BoundaryRules"]}];
        system = Reduce[system, Reals] // Quiet;(*TODO Quiet. We use Quiet because of: Reduce::ratnz: Reduce was unable to solve the system with inexact coefficients. The answer was obtained by solving a corresponding exact system and numericizing the result.*)
        If[system === True ||
        	Head[system] === Equal || 
        	(Head[system] === And && DeleteDuplicates[Head /@ system] === Equal),
        	(*Print["FixedStepX: Good job! \n"];*)
        	sol = Solve[system, Reals];
        	If[Length[sol] == 1,
        		sol = First @ sol,
        		Print["Not just one solution: ", sol];
        		sol = {};
        	],
        	Print["Possible multiple solutions \n", {system, rules}];
        	sol = {}; 
        ];
        Append[rules, sol] /. sol
    ]
EqEliminatorX3[{AA_, rules_List}] :=
    Module[ {EE, AA2, newrules, Asrules = Association @ rules, number},
    	(*Print["rules List"];*)
    	(*Print["System Heads: ",List @@ DeleteDuplicates@(Head[#] & /@ AA)];*)
    	EE = Select[AA, Head[#] === Equal &];
        AA2 = Select[AA, Head[#] =!= Equal &];
        (*Print["Equalities: ", EE];*)
        newrules = Solve[EE, Reals]//Quiet;
        number = Length @ newrules;
        If[number > 1,
        	Print["EqEliminatorX2: Number of solutions: ", Length @ newrules]];
        AssociateTo[Asrules, newrules // First];
        		(*TODO try without booleanconvert and reduce here*)        
        {BooleanConvert[Reduce[AA2/.Asrules, Reals] // Quiet,"CNF"], Asrules/. Asrules} (*Quiet : Reduce::ratnz: Reduce was unable to solve the system with inexact coefficients. The answer was obtained by solving a corresponding exact system and numericizing the result.*) 
    ]
    
EqEliminatorX3[{AA_, rules_Association}] := 
Module[{EE, AA2 = AA, newrules, Asrules = rules, number},
    	(*Print["EqEliminartorX3: Count! "];*)
        Which[AA2 === True, 
    		{AA2, Asrules},
        	AA2 === False,
        	Print["The system is False. Returning the inputs."];
        	{AA2, Asrules},
        	True,
        	EE = Select[AA2, Head[#] === Equal &];(*TODO what if there is no equality in AA2?*)
        	If[EE === True,
        		{AA2, Asrules},
        		newrules = Solve[EE, Reals]//Quiet;
        		(*Print[newrules];*)
        		Which[number > 1,
        			Print["EqEliminatorX3: Number of solutions: ", Length @ newrules, "\nThe system is: ", AA2],
        			number < 1,
        			Print["No solution for the system: ", EE];
        		];
        		AssociateTo[Asrules, newrules // First // Quiet];
        		(*Print[Asrules];*)
        			{BooleanConvert[Reduce[AA2/.Asrules, Reals] // Quiet,"CNF"], Asrules/. Asrules} (*TODO return message and STOP if Total @ Values[Asrules] is not numeric*)
        	]
        ]
    ]
 
 FixedSolverStepX3[Eqs_Association][rules_] :=
    Module[ {system, nonlinear, newsolve},
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[(Eqs["Nrhs"] /. rules),10^-6]}]);
        system = Eqs["EqAllAll"] && nonlinear;
        (*give system with sol to structural equations to equalities eliminator*)
        newsolve = FixedPoint[EqEliminatorX3, {system, Eqs["BoundaryRules"]}];
(*        Print["FixedX3: newsolve: ", {newsolve[[1]], Values@KeySort@newsolve[[2]]}];*)
        Print["FixedX3: next step: ", And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[(Eqs["Nrhs"] /. newsolve[[2]]),10^-6]}])];
        (*Print[N@Total @ Chop[(Eqs["Nrhs"] /. newsolve[[2]]),10^-6]];*)
        If[Head @ N @ Total @ Chop[(Eqs["Nrhs"] /. newsolve[[2]]),10^-6] === Real, 
        newsolve[[2]], Throw[rules ]//Quiet]
    ]


End[]
(*
EqEliminatorX2[{AA_, rules_Association}] := 
    Module[ {EE, AA2, newrules, Asrules = rules, number},
        EE = Select[AA, Head[#] === Equal &];
        If[EE === True,
        	AA2 = BooleanConvert[Reduce[AA, Reals] // Quiet,"CNF"], (*Quiet : Reduce::ratnz: Reduce was unable to solve the system with inexact coefficients. The answer was obtained by solving a corresponding exact system and numericizing the result.*)
        	AA2 = Select[AA, Head[#] =!= Equal &]
        ];
        	newrules = Solve[EE, Reals]//Quiet;
        	number = Length @ newrules;
        	If[number > 1,
        		Print["EqEliminatorX2: Number of solutions: ", Length @ newrules],
        		AssociateTo[Asrules, newrules // First];
        	{AA2, Asrules}/. Asrules 
        	]
        ];
        (*Print["System Heads: ",List @@ DeleteDuplicates@(Head[#] & /@ AA2)];*)
    	If[AA2 === True, 
        	{AA2, Asrules},
        	AssociateTo[Asrules, Solve[EE, Reals] // First // Quiet];
        	{AA2, Asrules} /. Asrules
        ]
    ]

*)    