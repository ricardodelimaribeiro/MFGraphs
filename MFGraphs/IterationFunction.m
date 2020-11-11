(* ::Package:: *)

(* Wolfram Language package *)

FixedReduceX1::usage = 
"FixedReduceX1[MFGEquations][rules] gives the solution for the system when the RHS of the nonlinear equations have rules substituted in them.
	MFGEquations should have the Keys EqAllAll, BoundaryRules, Nlhs, Nrhs, TOL";

EqEliminatorX::usage = 
"EqEliminator[{system, rules}] extracts the equalities of the system, solves and appends them to rules, and returns a new pair {system, rules} /. rules. ";

ZAnd::usage =
""

Begin["`Private`"]
CleanAndReplace::usage =
"CleanAndReplace[{system,rules}] returns True if the system is approximately (difference in each equation is less than 10^(-10) ) solved by rules.";
CleanAndReplace[{system_,rules_}] :=
    (system /. Equal -> (zob[#1-#2]&)/. rules // Chop[#,10^(-12)]& ) /. zob -> (# == 0.&)
 
ZAnd[_, False] :=
    False;
ZAnd[False, _] :=
    False;
ZAnd[xp_, True] :=
    xp;
ZAnd[xp_, eq_Equal] :=
    With[ {sol = Solve[eq]},
        If[ sol === {},
            False,
            (xp /. First[sol]) && And @@ (First[sol] /. Rule -> Equal) // 
     Simplify
        ]
    ];
ZAnd[xp_, orxp_Or] :=
    (ZAnd[xp, #] & /@ orxp) // RemoveDuplicates;

ZAnd[xp_, andxp_And] :=
    With[ {fst = First[andxp], rst = Rest[andxp]},
        If[ Head[fst] === Equal, 
         
         (* first is an equality *)
            With[ {sol = Solve[fst]},
                If[ sol === {},
                    False,
                    ZAnd[(xp /. First[sol]) && And @@ (First[sol] /. Rule -> Equal) //
        Simplify, ReplaceSolution[ rst, First[sol]]]
                ]
            ], 
            (* first is an OR - there are no other alternatives *)
            With[ {
               sol1 = Solve[fst[[1]]], 
               sol2 = Solve[fst[[2]]]
               },
                If[ sol1 === False,
                    False,
                    ZAnd[(xp /. First[sol1]) && 
                       And @@ (First[sol1] /. Rule -> Equal) // Simplify, 
                     ReplaceSolution[ rst, First[sol1]]]
                ]
                 ||
                 If[ sol2 === False,
                     False,
                     ZAnd[(xp /. First[sol2]) && 
                        And @@ (First[sol2] /. Rule -> Equal) // Simplify, 
                      ReplaceSolution[ rst, First[sol2]]]
                 ]
            ] // RemoveDuplicates
        ]
    ]

ReplaceSolution[True, sol_] :=
    True;
ReplaceSolution[False, sol_] :=
    False;
ReplaceSolution[rst_And, sol_] :=
    And @@ (Simplify /@ ((List @@ rst) /. sol))
ReplaceSolution[rst_, sol_] :=
    rst /. sol


RemoveDuplicates[xp_And] :=
    DeleteDuplicates[Sort[xp]];
RemoveDuplicates[xp_Or] :=
    DeleteDuplicates[Sort[xp]];
RemoveDuplicates[xp_] :=
    xp;


(*XPSort[xp_And] := Sort[xp];
XPSort[xp_] := xp;



NewReduce[True] :=
    True;


NewReduce[system_] :=
    Module[ {result},
        result = ZAnd[Select[system, !((Head[#] === Or)||(Head[#]===Equal))&],Select[system, ((Head[#] === Or)||(Head[#]===Equal))&]];
        Print["NewReduce: ", Head[result], result];
        If[Head[result] === Or,
        	result = RemoveDuplicates[
        		Reduce[Select[#, (Head[#] =!= Equal) &], Reals] &&
    			Select[#, (Head[#] === Equal) &] & /@ result]//Quiet,
        result = Reduce[#, Reals]& /@ result // Quiet
        ];(*Reduce was unable to solve the system with inexact coefficients. The answer was obtained by solving a corresponding exact system and numericizing the result.*)
        Reduce @ result
        (**)
	];
*)


NewReduce[True] :=
    True;


NewReduce[system_] :=
    Module[ {result},
        result = ZAnd[Select[system, !((Head[#] === Or)||(Head[#]===Equal))&],Select[system, ((Head[#] === Or)||(Head[#]===Equal))&]]//DeleteDuplicates;
        result = Reduce[#, Reals]& /@ result // Quiet;(*Reduce was unable to solve the system with inexact coefficients. The answer was obtained by solving a corresponding exact system and numericizing the result.*)
        Simplify @ result
    ];

EqEliminatorX[{system_, rules_}] :=
    Module[ {EE, ON, newrules, rulesAss = Association[rules]},
        (*Print["EEX {system, rules}:", {system,rules}];*)
        Which[
            Head[system] === And,
                EE = Select[system, (Head[#] === Equal) &];
                ON = Select[system, (Head[#] =!= Equal) &];
                If[ EE === {},
                    {Reduce[system, Reals], rulesAss},
                    (*Print["EEX: ", EE];*)
                    newrules =  Solve[EE, Reals] // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)
                       (*Print["EEX: newrules before if: ", newrules];*)
                    If[ newrules === {},
                           (*Print["EEX: ", EE];*)
                        {ON && Simplify @ EE, rulesAss},(*{system, rulesAss},*)
                        newrules = First @ newrules;  
                        (*Print["EEX: {system, newrules}: ", {system, newrules}];*)
                        (*Print["EEX: Really solve? \n", EE, EE /. newrules];
                        Print["EEX: Clean and Replace: \n", CleanAndReplace[{system, newrules}]];
                        Print["EEX: Clean and Replace just equalities: \n", CleanAndReplace[{EE, newrules}]];
                        Print["EEX: just substitute and simplify: \n", Simplify @ (system /. newrules)];*)
                        If[ Simplify @ (EE /. newrules) === False,
                            {ON && CleanAndReplace[{EE, newrules}] , AssociateTo[rulesAss, newrules]} /. newrules,
                            {system, AssociateTo[rulesAss, newrules]} /. newrules
                        ]
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
                {Reduce[system, Reals] // Quiet, rulesAss}
        ]
    ];

FixedReduceX1[Eqs_Association][rules_] :=
    Module[ {nonlinear, system = Eqs["EqAllAll"], auxsys, auxsol},
        (*These two steps can be done in DataToEquations*)
        {auxsys, auxsol} = FixedPoint[EqEliminatorX, {system, Eqs["BoundaryRules"]}]; 
            (*Print["FRX1: Structural stuff: \n", {auxsys, auxsol//KeySort}];*)
        auxsys = NewReduce[auxsys];
            (*Print["FRX1: Structural stuff -> NewReduce: \n", auxsys];*)
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[#,Eqs["TOL"]]&/@(Eqs["Nrhs"]/. auxsol /. rules)}]);
            (*Print["FRX1: nonlinear: \n", nonlinear];*)
        auxsys = (auxsys && nonlinear) /. auxsol;  
            (*Print["FRX1: new nonlinear: \n", {auxsys, auxsol, auxsys/.auxsol}];*)  
            (*Print["FRX1: ", {auxsys /. Solve[auxsys,Reals], Solve[auxsys,Reals]}];*)
        {auxsys, auxsol} = FixedPoint[EqEliminatorX, {auxsys, auxsol}];
        (*Print["FRX1: Structural stuff + nonlinear sol: \n", auxsol];*)
            (*Print["FRX1: Structural stuff + nonlinear: \n", {auxsys, auxsol//KeySort}];*)
            (*Print["FRX1: Structural stuff + nonlinear -> Reduce: \n", auxsys, Reduce[auxsys,Reals]];*)
        If[ auxsys === False,
            (*Print["rules :", rules];*)
            Print["FRX1: The last system in the iteration was inconsistent.\nThrowing the last feasible solution."];
            (*auxsol = rules;*)
            Throw[rules],
    (*        Print["FRX1: The relative Error on the nonlinear terms is ", Norm[(Eqs["Nlhs"] - Eqs["Nrhs"])/Eqs["Nlhs"] /. auxsol]]//Quiet;
            Print["FRX1: The relative Error on the nonlinear terms is ", Norm[(Eqs["Nlhs"] - Eqs["Nrhs"])/Eqs["Nrhs"] /. auxsol]]//Quiet;*)
            Print["FRX1: The Minimum between the Error and the relative Errors on the nonlinear terms is ", 
                Min[
                	Select[{
                		Norm[(Eqs["Nlhs"] - Eqs["Nrhs"]) /. auxsol],
                		(Norm[Eqs["Nlhs"] - Eqs["Nrhs"]]/Norm[Eqs["Nlhs"]]) /. auxsol,
                		(Norm[Eqs["Nlhs"] - Eqs["Nrhs"]]/Norm[Eqs["Nrhs"]]) /. auxsol
                		(*Norm[(Eqs["Nlhs"] - Eqs["Nrhs"])/Eqs["Nrhs"] /. auxsol]*)}, !(# === Indeterminate) &]]]//Quiet;
        ];
        Simplify /@ auxsol
    ];


End[]


(*FixedReduceX1[Eqs_Association][rules_] :=
    Module[ {nonlinear, system = Eqs["EqAllAll"], auxsys, auxsol},
    	(*Print["FRX1: first EEX..."];*)
    	{auxsys, auxsol} = FixedPoint[EqEliminatorX, {system, Eqs["BoundaryRules"]}];
        (*Print["FRX1: Structural stuff: \n", {auxsys, auxsol//KeySort}];*)
        auxsys = Reduce[auxsys,Reals];
        (*Print["FRX1: Structural stuff -> Reduce: \n", auxsys];*)
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[#,Eqs["TOL"]]&/@(Eqs["Nrhs"] /. rules)}]);
        (*Print["FRX1: nonlinear: \n", nonlinear];*)
        auxsys = auxsys && nonlinear /. auxsol;        
        (*Print["FRX1: second EEX..."];*)
        {auxsys, auxsol} = FixedPoint[EqEliminatorX, {auxsys, auxsol}];
        (*Print["FRX1: Structural stuff + nonlinear: \n", {auxsys, auxsol//KeySort}];*)
        (*Print["FRX1: Structural stuff + nonlinear -> Reduce: \n", Reduce[auxsys,Reals]];*)
        Print["FRX1: Error on the nonlinear terms: ", Norm[Eqs["Nlhs"] - Eqs["Nrhs"] /. auxsol]];
        auxsol
        (*Reduce[Eqs["AllAll"] && nonlinear, Reals]*)
    ];
*)


(*

  
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

    
*)