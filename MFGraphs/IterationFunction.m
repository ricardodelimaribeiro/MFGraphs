(* ::Package:: *)

(* Wolfram Language package *)

FixedReduceX1::usage = 
"FixedReduceX1[MFGEquations][rules] gives the solution for the system when the RHS of the nonlinear equations have rules substituted in them.
	MFGEquations should have the Keys EqAllAll, BoundaryRules, Nlhs, Nrhs, TOL";

CleanEqualities::usage = 
"CleanEqualities[{system,rules}] returns an equivalent pair of system and rules. The equalities are solved and substituted into the system and added to rules."

EqEliminatorX::usage = 
"EqEliminatorX[{system, rules}] extracts the equalities of the system, solves and appends them to rules, and returns a new pair {system, rules} /. rules. ";

ZAnd::usage =
""

NewReduce::usage = 
""

Begin["`Private`"]
CleanAndReplace::usage =
"CleanAndReplace[{system,rules}] returns True if the system is approximately (difference in each equation is less than 10^(-12) ) solved by rules.";
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
		(xp /. First[sol]) && And @@ (First[sol] /. Rule -> Equal) // Simplify
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
        			ZAnd[(xp /. First[sol]) && And @@ (First[sol] /. Rule -> Equal) // Simplify, ReplaceSolution[ rst, First[sol]]]
            ]
        ], 
        (* first is an OR - there are no other alternatives *)
        With[ {sol1 = Solve[fst[[1]]], sol2 = Solve[fst[[2]]]},
        		If[ sol1 === False,
        			False,
        			ZAnd[(xp /. First[sol1]) && And @@ (First[sol1] /. Rule -> Equal) // Simplify, ReplaceSolution[ rst, First[sol1]]]
        		]
            ||
            If[ sol2 === False,
            		False,
            		ZAnd[(xp /. First[sol2]) && And @@ (First[sol2] /. Rule -> Equal) // Simplify, ReplaceSolution[ rst, First[sol2]]]
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


NewReduce[True] :=
True;

NewReduce[False] :=
False;

NewReduce[system_] :=
Module[ {result},
	result = ZAnd[Select[system, !((Head[#] === Or)||(Head[#]===Equal))&],Select[system, ((Head[#] === Or)||(Head[#]===Equal))&]]//DeleteDuplicates;
	result = Reduce[#, Reals]& /@ result // Quiet;(*Reduce was unable to solve the system with inexact coefficients. The answer was obtained by solving a corresponding exact system and numericizing the result.*)
	Simplify @ result
]

EqEliminatorX[{system_, rules_}] :=
Module[ {EE, ON, newrules, rulesAss = Association[rules]},
	Which[
		Head[system] === And, 
		(*separete equalities from the rest*)
		EE = Select[system, (Head[#] === Equal) &];
		ON = Select[system, (Head[#] =!= Equal) &];
		If[ EE === {}, 
			(*there were no equalities: just reduce the system*)
			{Reduce[system, Reals], rulesAss},
			(*solve the equalities*)
			newrules =  Solve[EE, Reals] // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)(*TODO include variables to solve: this way we leave the switching costs unsolved.*)
			If[ newrules === {},
				(*equalities have no solution: simplify them*)
				{ON && Simplify @ EE, rulesAss},
                (*select the first set of solutions*)
                newrules = First @ newrules;  (*TODO: do we have more solutions?*)
                If[ Simplify @ (EE /. newrules) === False,
                    	(*maybe there were numeric errors!*)
                    {ON && CleanAndReplace[{EE, newrules}] , AssociateTo[rulesAss, newrules]} /. newrules,
                    (*solution is exact!*)
                    {system/. newrules, Simplify/@(AssociateTo[rulesAss, newrules]/. newrules)} 
                ]
            ]
        ],
        Head[system] === Equal,
        (*one last equality: this is rare! solutions are always unique and satisfy the equation*)
        newrules = First @ Solve[system, Reals];
        (*Simplify carries on computations like 1+1 in the Association*)
        {system/. newrules, Simplify/@(AssociateTo[rulesAss, newrules]/. newrules)},
        system === True,
        (*we need this because this function is used with FixedPoint. 
        Makes sure rules are simplified.*)
        {system, Simplify @ rulesAss},
        True,
        (*Maybe the Head is Or or NonNegative (could it be something else?)... 
        Just reduce the system while simplifyung the rules.*)
        {Reduce[system, Reals] // Quiet, Simplify @ rulesAss}
	]
];

CleanEqualities[{system_,rules_}] := FixedPoint[EqEliminatorX, {system,rules}];

FixedReduceX1[Eqs_Association][rules_] :=
Module[ {nonlinear, system = Eqs["nocasesystem"], auxsys, auxsol},
	alpha = Lookup[Eqs["Data"], "alpha", Function[{edge}, 1](*critical congestion is the default.*)];
	(*These two steps can be done in DataToEquations, but the NewReduce was somewhat slow...*)
    {auxsys, auxsol} = {system, Eqs["nocaserules"]}; 
    (*Print["FRX1: Structural stuff: \n", {auxsys, auxsol//KeySort}];*)
    (*auxsys = NewReduce[auxsys]; Done in DataToEquations*)
    (*Print["FRX1: Structural stuff -> NewReduce: \n", auxsys];*)
    (*Print["FRX1: edge -> hamiltonian: ", Eqs["Nrhs"]];*)
    nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[#,Eqs["TOL"]]&/@(Eqs["Nrhs"]/. auxsol /. rules)}]);
    (*Print["FRX1: nonlinear: \n", nonlinear];*)
    auxsys = (auxsys && nonlinear) /. auxsol;  
    (*Print["FRX1: new nonlinear: \n", {auxsys, auxsol, auxsys/.auxsol}];*)  
    (*Print["FRX1: ", {auxsys /. Solve[auxsys,Reals], Solve[auxsys,Reals]}];*)
    {auxsys, auxsol} = CleanEqualities[{auxsys, auxsol}];
    (*Print["FRX1: Structural stuff + nonlinear sol: \n", auxsol];*)
    (*Print["FRX1: Structural stuff + nonlinear: \n", {auxsys, auxsol//KeySort}];*)
    (*Print["FRX1: Structural stuff + nonlinear -> Reduce: \n", auxsys, Reduce[auxsys,Reals]];*)
    If[ auxsys === False,
    		(*precision issues...?but it could be the case that the "real" solution is repulsive/unstable.*)
    		(*Print["rules :", rules];*)
    		Print["FRX1: The last system in the iteration was inconsistent.\nReturning the last feasible solution.\nConsider that the solution may be unstable."];
    		(*auxsol = rules;*)
    		Return[rules],
    		(*Print["FRX1: The relative Error on the nonlinear terms is ", Norm[(Eqs["Nlhs"] - Eqs["Nrhs"])/Eqs["Nlhs"] /. auxsol]]//Quiet;
    		Print["FRX1: The relative Error on the nonlinear terms is ", Norm[(Eqs["Nlhs"] - Eqs["Nrhs"])/Eqs["Nrhs"] /. auxsol]]//Quiet;*)
    		Print["FRX1: The Minimum between the Error and the relative Errors on the nonlinear terms is ", 
    			Min[Select[{Norm[(Eqs["Nlhs"] - Eqs["Nrhs"]) /. auxsol],
    				(Norm[(Eqs["Nlhs"] - Eqs["Nrhs"])/Norm[Eqs["Nlhs"]]]) /. auxsol,
    				(Norm[(Eqs["Nlhs"] - Eqs["Nrhs"])/Norm[Eqs["Nrhs"]]]) /. auxsol
    				(*Norm[(Eqs["Nlhs"] - Eqs["Nrhs"])/Eqs["Nrhs"] /. auxsol]*)}, !(# === Indeterminate) &]
    			]
    		]//Quiet;
    ];
    Simplify /@ auxsol
];
End[]