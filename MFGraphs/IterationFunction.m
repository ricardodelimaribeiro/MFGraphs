(* ::Package:: *)

(* Wolfram Language package *)

Solver::usage =
"Solver[Eqs_Association][alpha] solves the network problem!"

FixedReduceX1::usage = 
"FixedReduceX1[MFGEquations][rules] gives the solution for the system when the RHS of the nonlinear equations have rules substituted in them.
	MFGEquations should have the Keys EqAllAll, BoundaryRules, Nlhs, Nrhs, TOL";

CleanEqualities::usage = 
"CleanEqualities[{system,rules}] returns an equivalent pair of system and rules without \"loose\" equalities. The equalities are solved and substituted into the system and added to rules."

FixedReduceN::usage = 
""

EqEliminatorX::usage = 
"EqEliminatorX[{system, rules}] extracts the equalities of the system, solves and appends them to rules, and returns a new pair {system, rules} /. rules. ";

EqEliminator::usage = 
"EqEliminator[{system, rules}] extracts the equalities of the system, solves and appends them to rules, and returns a new pair {system, rules} /. rules. ";

ZAnd::usage =
""

NewReduce::usage = 
""

ReZAnd::usage =
""

Begin["`Private`"]
CleanAndReplace::usage =
"CleanAndReplace[{system,rules}] returns True if the system is approximately (difference in each equation is less than 10^(-12) ) solved by rules.";
CleanAndReplace[{system_,rules_}] :=
    (system /. Equal -> (zob[#1-#2]&)/. rules // Chop(*[#,10^(-12)]&*) ) /. zob -> (# == 0.&)
 
ZAnd[_, False] :=
    False;

ZAnd[xp_, leq_LessEqual] :=
    (
        Simplify[xp && leq]
    )

ZAnd[xp_, geq_GreaterEqual] :=
    (    
        Simplify[xp && geq]
    )


ZAnd[xp_, ineq_Inequality] :=
    Simplify[xp && ineq]

ZAnd[False, _] :=
    False

ZAnd[xp_, True] :=
    xp

ZAnd[xp_, eq_Equal] :=
    With[ {sol = Solve[eq]},
        If[ sol === {},
            False,
            (xp /. First[sol]) && And @@ (First[sol] /. Rule -> Equal) // Simplify
        ]
    ]

ZAnd[xp_, orxp_Or] :=
    (    
        (*Print["ZAnd OR:\n",orxp];*)
        (ZAnd[xp, #] & /@ orxp) // RemoveDuplicates
    )
ZAnd[xp_, andxp_And] :=
    With[ {fst = First[andxp], rst = Rest[andxp]},
        Which[
            Head[fst] === Or,
            ReZAnd[xp,rst,#]& /@ fst (*// RemoveDuplicates*),
            True,
            ReZAnd[xp, rst, fst]
        ]
    ]
    
ReZAnd[xp_,rst_,fst_Equal] :=
    Module[ {fsol = First @ Solve @ fst},
        ZAnd[(xp /. fsol) && And @@ (fsol /. Rule -> Equal) // Simplify, ReplaceSolution[rst, fsol]]
    ]
  
ReZAnd[xp_,rst_,fst_Inequality] :=
	ZAnd[xp && fst, rst]
	
ReZAnd[xp_,rst_,fst_LessEqual] :=
    ZAnd[xp && fst, rst]

ReZAnd[xp_,rst_,fst_GreaterEqual] :=
    ZAnd[xp && fst, rst]

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

NewReduce[s_Or] :=
    s

NewReduce[True] :=
    True;

NewReduce[False] :=
    False;
NewReduce[s_Inequality] :=
    s
NewReduce[system_] :=
    Module[ {result},
        result = ZAnd[Select[system, !((Head[#] === Or)||(Head[#]===Equal))&],Select[system, ((Head[#] === Or)||(Head[#]===Equal))&]];
        (*Simplify @ result*)
        If[ result === False,
            False,
            (*BooleanConvert[result//DeleteDuplicates,"CNF"]*)
            result
        ]
    ]

EqEliminatorX[{system_, rules_}] :=
    Module[ {EE, ON, newrules, rulesAss = Association[rules]},
        Which[
            Head[system] === And, 
            (*separete equalities from the rest*)
            EE = Select[system, (Head[#] === Equal) &];
            If[ EE === {}, 
                (*there were no equalities: *)
                {system, rulesAss},
                (*solve the equalities*)
                newrules =  Solve[EE] // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)(*TODO include variables to solve: this way we leave the switching costs unsolved.*)
                ON = Select[system, (Head[#] =!= Equal) &];
                If[ newrules === {},
                    Print["EqEliminatorX: No solutions: \n", EE];
                    (*equalities have no solution: simplify them*)
                    {ON && Simplify @ EE, rulesAss},
                    (*select the first set of solutions*)
                    newrules = First @ newrules;  (*TODO: do we have more solutions? No! Linear equations...*)
                    If[ Simplify @ (EE /. newrules) === False,
                            (*maybe there were numeric errors!*)
                        {ON && CleanAndReplace[{EE, newrules}] , AssociateTo[rulesAss, newrules]} /. newrules,
                        (*solution is exact!*)
                        {ON /. newrules, Simplify/@(AssociateTo[rulesAss, newrules]/. newrules)}
                    (*changed system for ON*)
                    ]
                ]
            ],
            Head[system] === Equal,
            (*one last equality: this is rare! solutions are always unique and satisfy the equation*)
            newrules = First @ Solve[system];
            (*Simplify carries on computations like 1+1 in the Association*)
            {system/. newrules, Simplify/@(AssociateTo[rulesAss, newrules]/. newrules)},
            system === True,
            (*we need this because this function is used with FixedPoint. 
            Makes sure rules are simplified.*)
            {system, Simplify @ rulesAss},
            True,
            (*Maybe the Head is Or or NonNegative (could it be something else?)... *)
            {system, Simplify @ rulesAss}
        ]
    ];

EqEliminator[{True, rules_ }] :=
    {True, Simplify /@ rules}

EqEliminator[{system_Or, rules_ }] :=
    ((*Print[system];
    Print[Simplify@system];*)
       {Simplify @ system, Simplify /@ rules})
       
EqEliminator[{system_GreaterEqual, rules_ }] :=
    ((*Print[system];
    Print[Simplify@system];*)
       {Simplify @ system, Simplify /@ rules})
       
EqEliminator[{system_LessEqual, rules_ }] :=
    ((*Print[system];
    Print[Simplify@system];*)
       {Simplify @ system, Simplify /@ rules})
       
EqEliminator[{system_Inequality, rules_ }] :=
    ((*Print[system];
    Print[Simplify@system];*)
       {Simplify @ system, Simplify /@ rules})

       
EqEliminator[{False,rules_}] :=
    {False, rules}


EqEliminator[{system_, rules_}] :=
    Module[ {EE, ON, newrules, rulesAss = Association[rules]},
        (*separete equalities from the rest*)
        EE = Select[system, (Head[#] === Equal) &];
        (*Print["EE:\n", EE];*)
        If[ EE === False,
            Print["EL: ", system];
            Return[{system, rulesAss}]
        ];
        newrules =  Solve[EE] // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)(*TODO include variables to solve: this way we leave the switching costs unsolved.*)
        ON = Select[system, (Head[#] =!= Equal) &];
        If[ Length[newrules] == 1,
            newrules = First @ newrules;
            {BooleanConvert[ON//.newrules,"CNF"], Simplify /@ (AssociateTo[rulesAss, newrules]//. newrules)}, (*changed system for ON*)
            Print["Rules for the equalities: \n", newrules,"\nSystem: \n", EE];
            {system,rulesAss}
        ]
    ];

(*CleanEqualities[{system_,rules_}] := FixedPoint[EqEliminatorX, {system,rules}];*)
CleanEqualities[{system_, rules_}] :=
    FixedPoint[EqEliminator, {system,rules}];

FixedReduceX1[Eqs_Association][rules_] :=
    Module[ {nonlinear, auxsys, auxsol, error, TOL = 10^-10},
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], RoundValues[Eqs["Nrhs"]/. rules]}]);
        (*{nonlinear, auxsol} = CleanEqualities[{nonlinear,{}}];*)
        auxsol = First[Solve[nonlinear]//Quiet];
        auxsys = (Eqs["EqAllAll"] && nonlinear) /. auxsol;
        {auxsys, auxsol} = CleanEqualities[{auxsys, auxsol}];
        (*auxsys = NewReduce[BooleanConvert[auxsys,"CNF"]];*)
        auxsys = NewReduce[auxsys];
        {auxsys, auxsol} = CleanEqualities[{auxsys, auxsol}];
        auxsol = Simplify @ auxsol;
        If[ auxsys === False,
            (*precision issues...?but it could be the case that the "real" solution is repulsive/unstable.*)
            Print["FRX1: The last system in the iteration was inconsistent.\nReturning the last feasible solution.\nConsider that the solution may be unstable."];
            Return[rules],
            Print["The error (1-Norm of LHS-RHS) is ", Norm[(Eqs["Nlhs"] - Eqs["Nrhs"]) /. auxsol, 1]];
        ];
        error = Norm[(Eqs["Nlhs"] - Eqs["Nrhs"]) /. auxsol, 1];
        If[ error<TOL,
            Print["The error is ", error//ScientificForm, " which is less than " , TOL//N//ScientificForm];
            Throw[auxsol]
        ];
        auxsol
    ];

Solver[Eqs_Association][alpha_] :=
    Module[ {rules0,FFR, timeFFR},
        If[ alpha == 1 || alpha == 1.,
            Module[ {js = Lookup[Eqs, "js", Return[]]},
                rules0 = AssociationThread[js,0& /@ js];
                {timeFFR, FFR} = AbsoluteTiming@Catch@FixedPoint[FixedReduceX1[Eqs], rules0, 1];
                Print["It took ", timeFFR, " seconds to solve!\n\tThe system is ",Simplify[Eqs["EqAllAll"] && Eqs["EqCriticalCase"]/.FFR]];
                FFR
            ],
            {timeFFR, FFR} = AbsoluteTiming[Catch@FixedPoint[FixedReduceX1[Eqs], rules0, 2] // N // KeySort];
            FFR
        ]
    ]

FixedReduceN[Eqs_Association][rules_] :=
    Module[ {nonlinear, system = Eqs["EqAllAll"], auxsys, auxsol},
        (*alpha = Lookup[Eqs["Data"], "alpha", Function[{edge}, 1](*critical congestion is the default.*)];*)
        Print[alpha];
        (*These two steps can be done in DataToEquations, but the NewReduce was somewhat slow...*)
        {auxsys, auxsol} = {system, <||>}; 
        (*Print["FRX1: Structural stuff: \n", {auxsys, auxsol//KeySort}];*)
        (*auxsys = NewReduce[auxsys]; Done in DataToEquations*)
        (*Print["FRX1: Structural stuff -> NewReduce: \n", auxsys];*)
        (*Print["FRX1: edge -> hamiltonian: ", Eqs["Nrhs"]];*)
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], Chop[#,Lookup[Eqs,"TOL",10^10]]&/@(Eqs["Nrhs"]/. auxsol /. rules)}]);
        {nonlinear,auxsol} = CleanEqualitiesN[{nonlinear,auxsol}];
        (*Print["FRX1: nonlinear: \n", auxsol];*)
        auxsys = (auxsys && nonlinear) /. auxsol;
        {auxsys, auxsol} = CleanEqualitiesN[{auxsys, auxsol}];
        (*Print["FRX1: Structural stuff + nonlinear sol before NewReduce: \n", {auxsys,auxsol}];*)
        auxsys = NewReduce[BooleanConvert[auxsys,"CNF"]];
        {auxsys, auxsol} = CleanEqualitiesN[{auxsys, auxsol}];
        auxsol = Simplify @ auxsol;
        (*Print["FRX1: Structural stuff + nonlinear sol: \n", {auxsys,auxsol}];*)
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
            Print[(Eqs["Nlhs"] - Eqs["Nrhs"]) /. auxsol];
        ];
        auxsol
    ];

End[]