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


NewReduce[True] :=
    True;
NewReduce[False] :=
False;

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
                        	(*Print["in if, true"];*)
                            {ON && CleanAndReplace[{EE, newrules}] , AssociateTo[rulesAss, newrules]} /. newrules,
                            (*Print["in if, false"];*)
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