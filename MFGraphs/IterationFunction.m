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

FixedReduceX1::usage = 
"";

EqEliminatorX::usage = 
"";

EqEliminatorX2::usage = 
"Improving the Eliminator.";

EqEliminatorX3::usage = 
"Improving the Eliminator.";

Begin["`Private`"]

   (*TODO critical congestion solution tester: replace "solution" on the equations. Look at the output to figure out if it works!! *)
   (*in the nonlinear problem, look for approximate equalities.*)
  
EqEliminatorX2[{AA_, rules_List}] :=
    Module[ {EE, AA2, newrules, Asrules = Association @ rules, number},
    	EE = Select[AA, Head[#] === Equal &];
        AA2 = Select[AA, Head[#] =!= Equal &];
        newrules = Solve[EE, Reals]//Quiet;
        number = Length @ newrules;
        If[number > 1,
        	Print["EqEliminatorX2: Number of solutions: ", Length @ newrules]];
        AssociateTo[Asrules, newrules // First];
        {AA2, Asrules}/. Asrules 
    ]
    
EqEliminatorX2[{AA_, rules_Association}] := 
Module[{EE, AA2, newrules, Asrules = rules, number},
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
        		number = Length[newrules];
        		Which[number > 1,
        			Print["EqEliminatorX2: Number of solutions: ", Length @ newrules, "\nThe system is: ", AA2],
        			number < 1,
        			Print["No solution for the system: ", EE];
        		];
        		AssociateTo[Asrules, newrules // First // Quiet];
        		{AA2, Asrules} /. Asrules (*TODO return message and STOP if Total @ Values[Asrules] is not numeric*)
        	]
        ]
    ] 

FixedSolverStepX2[Eqs_Association][rules_] :=
    Module[ {system, nonlinear, newsolve},
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[(Eqs["Nrhs"] /. rules),10^-6]}]);
        system = Eqs["EqAllAll"] && nonlinear;
        newsolve = FixedPoint[EqEliminatorX2, {system, Eqs["BoundaryRules"]}];
        Print["FixedX2: next step: ", And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[(Eqs["Nrhs"] /. newsolve[[2]]),10^-6]}])];
        newsolve[[2]]
    ]
    
EqEliminatorX[{system_, rules_}] := (*if it solves, system becomes True. rules is a list of rules or an association.*)
Module[{EE, newrules, rulesAss = Association[rules]},
    Which[
    	Head[system] === And,
    		EE = Select[system, (Head[#] === Equal) &];
            If[EE === {},
            	{Reduce[system, Reals], rulesAss},
               		newrules =  Solve[EE, Reals] // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)
               		If[newrules === {},
               			{system, rulesAss},
               			newrules = First @ newrules;  
               			{system, AssociateTo[rulesAss, newrules]}/.newrules
               		]
            ],
        Head[system] === Equal,
            newrules = Solve[system, Reals];
            Print["EEX: newrules when Head is Equal: ", newrules];
            newrules = First @ newrules;
            {system, AssociateTo[rulesAss, newrules]} /. newrules,        
        system === True,
        	{system, rulesAss},
        True,
        	(*Print["EEX: System is not True and Head of system is NOT And nor Equal. Returning the input:\n", {system, rulesAss}];
        	Print["EEX: Testing Reduce:\n", {Reduce[system, Reals], rulesAss}];*)
        	{Reduce[system, Reals] // Quiet, rulesAss}
    ]
];

FixedReduceX1[Eqs_Association][rules_] :=
    Module[ {nonlinear, system = Eqs["EqAllAll"], auxsys, auxsol},
    	(*Print["FRX1: first EEX..."];*)
    	{auxsys, auxsol} = FixedPoint[EqEliminatorX, {system, Eqs["BoundaryRules"]}];
        (*Print["FRX1: Structural stuff: \n", {auxsys, auxsol//KeySort}];*)
        auxsys = Reduce[auxsys,Reals];
        (*Print["FRX1: Structural stuff -> Reduce: \n", auxsys];*)
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[#,10^-3]&/@(Eqs["Nrhs"] /. rules)}]);
        (*Print["FRX1: nonlinear: \n", nonlinear];*)
        auxsys = auxsys && nonlinear /. auxsol;        
        (*Print["FRX1: second EEX..."];*)
        {auxsys, auxsol} = FixedPoint[EqEliminatorX, {auxsys, auxsol}];
        (*Print["FRX1: Structural stuff + nonlinear: \n", {auxsys, auxsol//KeySort}];*)
        (*Print["FRX1: Structural stuff + nonlinear -> Reduce: \n", Reduce[auxsys,Reals]];*)
        Print[Norm[Eqs["Nlhs"] - Eqs["Nrhs"] /. auxsol]];
        auxsol
        (*Reduce[Eqs["AllAll"] && nonlinear, Reals]*)
    ];

FixedSolverStepX1[Eqs_Association][rules_] :=
    Module[ {system = Eqs["EqAllAll"], nonlinear, sol, nsol},
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[(Eqs["Nrhs"] /. rules),10^-6]}]);
        nsol = Solve[nonlinear, Reals] // Quiet // First // Association;

        (*begin: for debuggin purposes*)
        (*end*)

        (*Print["FSSX1: ", nonlinear];*)
        (*Print["FSSX1: nonlinear rules: ", nsol];*)
        system = system /. nsol;
        AssociateTo[nsol, Eqs["BoundaryRules"]] /. nsol;
        (*Print["FSSX1: initial rules for EqEliminatorX: ", nsol];*)
        {system, nsol} = FixedPoint[EqEliminatorX, {system,  nsol}];
        Print["FSSX1: system: ", system];
        system = Reduce[system, Reals] // Quiet;(*TODO Quiet. We use Quiet because of: Reduce::ratnz: Reduce was unable to solve the system with inexact coefficients. The answer was obtained by solving a corresponding exact system and numericizing the result.*)
        Print["FSSX1: reduced system: ", system];
        If[system === True ||
        	Head[system] === Equal || 
        	(Head[system] === And && DeleteDuplicates[Head /@ system] === Equal),
        	Print["FSSX1: If is true"];
        	sol = Solve[system, Reals] // Quiet;
        	If[Length[sol] == 1,
        		sol = First @ sol,
        		Print["FixedSolverStepX1: Not just one solution: ", sol];
        		sol = {};
        	],
        	(*Print["FSSX1: If is false"];*)
        	Print["FixedSolverStepX1: System is ", system, "\nReturning the input rules/association: \n", rules];
        	sol = {}; 
        ];
        Print["FSSX1: End of code: {system, nonlinear}\n", {system,nonlinear},
        	"\nFSSX1: End of code: {system, nonlinear} /. nsol\n",{system,nonlinear} /. nsol,
        	"\nFSSX1: End of code: {system, nonlinear} /. (AssociateTo[nsol, sol] /. sol)\n",{system,nonlinear} /. (AssociateTo[nsol, sol] /. sol)];   
        AssociateTo[nsol, sol] /. sol
    ]
EqEliminatorX3[{AA_, rules_List}] :=
    Module[ {EE, AA2, newrules, Asrules = Association @ rules, number},
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


(*EqEliminatorX[{system_, rules_}] := (*if it solves, system becomes True. rules is a list of rules or an association.*)
    Which[
    	Head[system] === And, 
        	Module[ {EE, newrules, rulesAss = Association[rules]},
        		EE = Select[system, (Head[#] === Equal) &];
            	If[EE === {}, 
            		{system, rulesAss},
               		newrules =  Solve[EE, Reals] // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)
               		If[newrules === {},
               			{system, rulesAss},
               			newrules = First @ newrules;  
               			{system, AssociateTo[rulesAss, newrules]}/.newrules
               		] 
            	]
        	],
        Head[system] === Equal,
        	Module[ {newrules, rulesAss = Association[rules]},
            	newrules = Solve[system, Reals];
            	Print[newrules];
            	newrules = First @ newrules;
            	{system, AssociateTo[rulesAss, newrules]}/.newrules
        	],        
        system === True,
        	{system, Association @ rules},
        True,
        	Print["Eliminator: System is not True and Head of system is NOT And or Equal. Returning the input."];
        	{system, Association @ rules}
    ];
*)