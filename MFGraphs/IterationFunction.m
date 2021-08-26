(* ::Package:: *)

(* Wolfram Language package *)

CriticalCongestionSolver::usage = 
"CriticalCongestionSolver[eq_Association] returns the critical congestion solution."

CriticalCongestionStep::usage =
""

Solver::usage =
"Solver[Eqs_Association][alpha] solves the network problem with the given alpha."

FixedReduceX1::usage = 
"FixedReduceX1[MFGEquations][rules] gives the solution for the system when the RHS of the nonlinear equations have rules substituted in them.
	MFGEquations should have the Keys EqAllAll, BoundaryRules, Nlhs, Nrhs, TOL";(*TODO review the keys...*)

CleanEqualities::usage = 
"CleanEqualities[{system,rules}] returns an equivalent pair of system and rules without \"loose\" equalities. 
The equalities are solved and substituted into the system and added to rules."

EqEliminator::usage = 
"EqEliminator[{system, rules}] extracts the equalities of the system, solves and appends them to rules, and returns a new pair {system, rules} /. rules. ";

ZAnd::usage =
"ZAnd[system, orsys] is a recursive function to reduce system&&orsys by checking the feasibility of the alternatives in orsys."

NewReduce::usage = 
"NewReduce[system] reduces the system using ZAnd"

ReZAnd::usage =
"ReZAnd[] "

RemoveDuplicates::usage =
"RemoveDuplicates[exp] Sort and DeleteDuplicates"

EliminateVarsStep::usage = 
""

EliminateVars::usage = 
""

Begin["`Private`"]
CleanAndReplace::usage =
"CleanAndReplace[{system,rules}] returns True if the system is approximately (difference in each equation is less than 10^(-10) ) solved by rules.";

CleanAndReplace[{system_,rules_}] :=
    (system /. Equal -> (zob[#1-#2]&)/. rules // Chop(*[#,10^(-12)]&*) ) /. zob -> (# == 0.&)
 
ZAnd[_, False] :=
    False

ZAnd[xp_, leq_LessEqual] :=
    Simplify[xp && leq]

ZAnd[xp_, geq_GreaterEqual] :=
    Simplify[xp && geq]

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
    (ZAnd[xp, #] & /@ orxp) // RemoveDuplicates
    
        (*TODO check that the structure is basically the same*)
ZAnd[xp_, andxp_And] :=
    With[ {fst = First[andxp], rst = Rest[andxp]},
        Which[
            Head[fst] === Or,
            ReZAnd[xp,rst] /@ fst (*// RemoveDuplicates*),
            True,
            ReZAnd[xp, rst, fst]
        ]
    ]
    
(*Operator form of ReZand*)    
ReZAnd[xp_,rst_] :=
    ReZAnd[xp,rst,#]&     
    
ReZAnd[xp_,rst_,fst_Equal] :=
    Module[ {fsol = First @ Solve @ fst // Quiet},
        ZAnd[(xp /. fsol) && And @@ (fsol /. Rule -> Equal) // Simplify, ReplaceSolution[rst, fsol]]
    ]

(*TODO fix this if needed*)
ReZAnd[xp_, rst_, fst_And] :=
    (
    Print["got And: ", fst]; 
    ReZAnd[xp,rst] /@ fst 
    )  
ReZAnd[xp_,rst_,fst_Inequality] :=
    ZAnd[xp && fst, rst]
    
ReZAnd[xp_,rst_,fst_LessEqual] :=
    ZAnd[xp && fst, rst]

ReZAnd[xp_,rst_,fst_GreaterEqual] :=
    ZAnd[xp && fst, rst]

ReplaceSolution[True, sol_] :=
    True

ReplaceSolution[False, sol_] :=
    False

ReplaceSolution[rst_And, sol_] :=
    Simplify /@ (rst /. sol)
    
ReplaceSolution[rst_, sol_] :=
    rst /. sol

RemoveDuplicates[xp_And] :=
    DeleteDuplicates[Sort[xp]];

RemoveDuplicates[xp_Or] :=
    DeleteDuplicates[Sort[xp]];
    
RemoveDuplicates[xp_] :=
    xp

NewReduce[s_Or] :=
    Simplify @ s

NewReduce[True] :=
    True

NewReduce[False] :=
    False
    
NewReduce[x_Equal] :=
    Simplify @ x
    
NewReduce[s_Inequality] :=
    s
    
NewReduce[x_GreaterEqual] :=
    x

NewReduce[x_LessEqual] :=
    x

NewReduce[x_Inequality] :=
    x

NewReduce[system_And] :=
    Module[ {result},
        result = ZAnd[Select[system, !((Head[#] === Or)||(Head[#]===Equal))&],Select[system, ((Head[#] === Or)||(Head[#]===Equal))&]];
        If[ result === False,
            False,
            result // DeleteDuplicates
        ]
    ]

EqEliminator[{True, rules_ }] :=
    {True, Simplify /@ rules}

EqEliminator[{system_Or, rules_ }] :=
    {Simplify @ system, Simplify /@ rules}
       
EqEliminator[{system_GreaterEqual, rules_ }] :=
    {Simplify @ system, Simplify /@ rules}
       
EqEliminator[{system_LessEqual, rules_ }] :=
    {Simplify @ system, Simplify /@ rules}
       
EqEliminator[{system_Inequality, rules_ }] :=
    {Simplify @ system, Simplify /@ rules}

EqEliminator[{False,rules_}] :=
    {False, rules}
    
EqEliminator[{system_, rules_List}] :=
    EqEliminator[{system, Association[rules]}]

EqEliminator[{system_Equal, rules_Association}] :=
Module[{newrules},
	newrules = First @ Solve[system];
	newrules = Association[newrules];
	{True, Join[rules/.newrules,newrules]}
]

EqEliminator[{system_And, rules_Association}] :=
    Module[ {EE, ON, newrules},
        (*separate equalities from the rest*)
        EE = Select[system, (Head[#] === Equal) &];

        If[ EE === True,
            Return[{system, rules}]
        ];

        newrules =  Solve[EE] // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)
        ON = Select[system, (Head[#] =!= Equal) &];
        newrules = First @ newrules;
        newrules = Join[rules, Association @ newrules];
        {ON /. newrules, newrules /. newrules}
    ]

CleanEqualities[system_] :=
    CleanEqualities[{system, {}}]

CleanEqualities[{system_, rules_}] :=
    FixedPoint[EqEliminator, {system,rules}]

FixedReduceX1[Eqs_Association][rules_] :=
    Module[ {nonlinear, auxsys, auxsol, error, TOL = 10^-10},
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], RoundValues[Eqs["Nrhs"]/. rules]}]);
        auxsol = First[Solve[nonlinear]//Quiet];
        auxsys = (Eqs["EqAllAll"] && nonlinear) /. auxsol;
        {auxsys, auxsol} = CleanEqualities[{auxsys, auxsol}];
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
    ]

Solver[Eqs_Association][alpha_] :=
    Module[ {rules0, FFR, timeFFR},
        If[ alpha == 1 || alpha == 1.,
            Module[ {js},
                js = Values@Lookup[Eqs, "jvars", Print["hey hey hey"];
                                                 Return[]];
                rules0 = AssociationThread[js,0& /@ js];
                {timeFFR, FFR} = AbsoluteTiming@Catch@FixedPoint[FixedReduceX1[Eqs], rules0, 1];
                Print["It took ", timeFFR, " seconds to solve!\n\tThe system is ",Simplify[Eqs["EqAllAll"] && Eqs["EqCriticalCase"]/.FFR]];
                FFR
            ],
            {timeFFR, FFR} = AbsoluteTiming[(*Catch@*)FixedPoint[FixedReduceX1[Eqs], CriticalCongestionSolver[Eqs][[2]], 5]];
            Print["It took ", timeFFR, " seconds to solve!\n\tThe system is ",Simplify[Eqs["EqAllAll"] /.FFR]];
            FFR
        ]
    ]

CriticalBundle[Data_Association] :=
    Module[ {d2e},
        d2e = D2E[Data];
        {d2e,CriticalCongestionSolver[d2e]}
    ]

CriticalCongestionSolver[D2E_Association] :=
    Module[ {system, rules, time,
    EqAllAll = Lookup[D2E, "EqAllAll", Print["No equations to solve."];
                                       Return[]], 
    EqCriticalCase = Lookup[D2E,"EqCriticalCase", Print["Critical case equations are missing."];
                                                  Return[]],
    InitRules = Lookup[D2E, "BoundaryRules", Print["Need boundary conditions"];],
    js =  Values @ Lookup[D2E, "jvars"],
    jts = Values @ Lookup[D2E, "jtvars"],
    us = Values @ Lookup[D2E, "uvars"]  
},
		Print["Number of variables: ",Length[Join[js,jts,us]]];
        {system,rules} = CleanEqualities[{EqCriticalCase && EqAllAll, InitRules}];   
        Print["Removed variables: ",Length[rules]];
        (*Lets solve for the u*)
        {time, {system, rules}} = AbsoluteTiming @ CleanEqualities[{(And @@ (NewReduce /@ (Select[system, Function[a, ! FreeQ[a, #]]] & /@ us)) // Simplify) && system, rules}];
        Print["Time to \"eliminate\" u variables: ",time," seconds and"];
        Print["Removed variables: ",Length[rules]];
        system = BooleanConvert[system,"CNF"];
        {system, rules} = CleanEqualities[{system, rules}];
        {system, rules} = FixedPoint[CriticalCongestionStep, {system, rules}, SameTest -> 
            (*(Equal[ #1[[2]], #2[[2]]]&)];*)
            (Equal[Length[ #1[[2]] ], Length[ #2[[2]] ] ] &) ];(*Try a different test by adding something like: || AllTrue[ (Function[x,FreeQ[x][ #2[[1]] ]] /@ js)]*)
	        If[ !(And @@ (FreeQ[#][system] & /@ js)),
            	Print["There are still currents in the system."];
            	system = Reduce[system, Reals];
            	Print[system];
            	{system, rules} = CleanEqualities[{system, rules}]
        	];
        rules = Simplify /@ rules /. rules;
        {system, rules} /. rules
    ]

CriticalCongestionStep[{sys_Or, rul_}] :=
    {BooleanConvert[sys,"CNF"], rul}

CriticalCongestionStep[{False, rul_}] :=
    {False, rul}

CriticalCongestionStep[{True, rul_}] :=
    {True, rul}

CriticalCongestionStep[{sys_,rul_}] :=
    (Print["Removed variables: ",Length[rul]];
     CleanEqualities[{NewReduce[BooleanConvert[sys,"CNF"](*sys*)], rul}])

EliminateVarsStep[{{system_, rules_}, {}}] :=
 {{system, rules}, {}}
 
EliminateVarsStep[{{system_, rules_}, vars_}] :=
 Module[{var, subsys, subsyscomplement, position , newsys, newrules},
  var = SelectFirst[vars, ! FreeQ[#][system] &];
  Print[var];
  If[Head[var] === Missing,
   {{system, rules}, {}},
   position = First@FirstPosition[vars, var];
   subsys = Select[system, Function[x, ! FreeQ[var][x]]];
   subsyscomplement = Select[system, Function[x, FreeQ[var][x]]];
   subsys = Simplify[subsys];
   newsys = subsys && subsyscomplement;
   {newsys, newrules} = CleanEqualities[{newsys, rules}];
   (*{{BooleanConvert[newsys,"CNF"],newrules},Drop[vars,
   position]}*)
   {{newsys, newrules}, Drop[vars, position]}
   ]
  ]
  
EliminateVars[{{system_, rules_}, vars_}] :=
 FixedPoint[EliminateVarsStep, {{system, rules}, vars}]
 
End[]