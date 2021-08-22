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
    True

ReplaceSolution[False, sol_] :=
    False

ReplaceSolution[rst_And, sol_] :=
    (*And @@ (Simplify /@ ((List @@ rst) /. sol))*)
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
    s

NewReduce[True] :=
    True;

NewReduce[False] :=
    False
    
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

EqEliminator[{system_, rules_}] :=
    Module[ {EE, ON, newrules, rulesAss = Association[rules]},
        (*separete equalities from the rest*)
        EE = Select[system, (Head[#] === Equal) &];
        (*Print["EE:\n", EE];*)
        If[ EE === False,
            Print["EL: ", system];
            Return[{system, rulesAss}]
        ];
        (*TODO include variables to solve: this way we leave the switching costs unsolved.*)
        newrules =  Solve[EE] // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)
        ON = Select[system, (Head[#] =!= Equal) &];
        If[ Length[newrules] == 1,
            newrules = First @ newrules;
            (*{BooleanConvert[ON//.newrules,"CNF"], Simplify /@ (AssociateTo[rulesAss, newrules]//. newrules)}, (*changed system for ON*)*)
            {ON//.newrules, Simplify /@ (AssociateTo[rulesAss, newrules]//. newrules)}, (*changed system for ON*)
            Print["Rules for the equalities: \n", newrules,"\nSystem: \n", EE];
            {system,rulesAss}
        ]
    ]

CleanEqualities[system_] :=
    CleanEqualities[{system, {}}]

CleanEqualities[{system_, rules_}] :=
    FixedPoint[EqEliminator, {system,rules}]

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
        CriticalCongestionSolver[d2e]
    ]

CriticalCongestionSolver[D2E_Association] :=
    Module[ {system, rules,
    EqAllAll = Lookup[D2E, "EqAllAll", Print["No equations to solve."];
                                       Return[]], 
    EqCriticalCase = Lookup[D2E,"EqCriticalCase", Print["Critical case equations are missing."];
                                                  Return[]],
    InitRules = Lookup[D2E, "BoundaryRules", Print["Need boundary conditions"];]  
},
		InitRules = Join[InitRules,First @ Solve[D2E["AllEq"]/.InitRules]];
        EqAllAll = EqAllAll /. InitRules;
        EqCriticalCase = EqCriticalCase /. InitRules;
        rules = First @ Solve @ EqCriticalCase//Quiet;
        system = EqAllAll/.rules;
        Print[system];
        {system,rules} = FixedPoint[CriticalCongestionStep, {system,Join[InitRules,rules]}];
        rules = Simplify /@ rules;
        Simplify /@ ({system, rules}/.rules)
    ]

CriticalCongestionStep[{sys_Or, rul_}] :=
    {Reduce[sys, Reals], rul}

CriticalCongestionStep[{False, rul_}] :=
    {False, rul}

CriticalCongestionStep[{sys_,rul_}] :=
    Module[ {system = sys,rules = rul},
        {system,rules} = CleanEqualities[{system, rules}];
        system = BooleanConvert[system,"CNF"];
        (*{system,rules} = CleanEqualities[{system, rules}];*)
        system = NewReduce[system];
        {system,rules}
    ]


End[]