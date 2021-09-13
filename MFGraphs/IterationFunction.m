(* ::Package:: *)

(* Wolfram Language package *)

CriticalCongestionSolver::usage = 
"CriticalCongestionSolver[eq_Association] returns the critical congestion solution."

CleanEqualitiesOperator::usage = 
"CleanEqualitiesOperator[D2E][{system,rules}] returns an equivalent pair of system and rules without \"loose\" equalities. 
The equalities are solved and substituted into the system and added to rules."

EqEliminator::usage = 
"EqEliminator[Eqs][{system, rules}] extracts the equalities of the system, solves and appends them to rules, and returns a new pair {system, rules} /. rules. ";

NewReduce::usage = 
"NewReduce[system] reduces the system using ZAnd"

Begin["`Private`"]
EliminateVarsStep::usage = 
""

EliminateVars::usage = 
""

ZAndEqual::usage =
"ZAndEqual[rules, xp] returns an association with the solution of the equalities"

NewCleanEqualities::usage =(*this strategy is too slow! maybe we can repair...*)
"NewCleanEqualities[{system, rules}] returns the non-equalities and the solution to the complementary linear system."

SimpleCrit::usage=
"SimpleCrit[{{syscrit, rules}, critus}] solves the critical congestion equations for the us."

RemoveDuplicates::usage =
"RemoveDuplicates[exp] Sort and DeleteDuplicates"

CriticalCongestionStep::usage =
""

ZAnd::usage =
"ZAnd[system, orsys] is a recursive function to reduce system&&orsys by checking the feasibility of the alternatives in orsys."

ReZAnd::usage =
"ReZAnd[] "


Solver::usage =
"Solver[Eqs_Association][alpha] solves the network problem with the given alpha."

FixedReduceX1::usage = 
"FixedReduceX1[MFGEquations][rules] gives the solution for the system when the RHS of the nonlinear equations have rules substituted in them.
	MFGEquations should have the Keys EqAllAll, BoundaryRules, Nlhs, Nrhs, TOL";(*TODO review the keys...*)


CleanAndReplace::usage =
"CleanAndReplace[{system,rules}] returns True if the system is approximately (difference in each equation is less than 10^(-10) ) solved by rules.";

CleanAndReplace[{system_,rules_}] :=
    (system /. Equal -> (zob[#1-#2]&)/. rules // Chop(*[#,10^(-12)]&*) ) /. zob -> (# == 0.&)
 
ZAnd[_, False] :=
    ((*Print["second is false"];*)
     False)

ZAnd[xp_, leq_LessEqual] :=
    Simplify[xp && leq]

ZAnd[xp_, geq_GreaterEqual] :=
    Simplify[xp && geq]

ZAnd[xp_, ineq_Inequality] :=
    Simplify[xp && ineq]

ZAnd[False, _] :=
    ((*Print["first is false"];*)
    
     False)

ZAnd[xp_, True] :=
    xp

ZAnd[xp_, eq_Equal] :=
    With[ {sol = Solve[eq]},
        If[ sol === {},
            False,
            (*Print[First@sol];*)
            (xp /. First[sol]) && And @@ (First[sol] /. Rule -> Equal) // Simplify
        ]
    ]

ZAnd[xp_, orxp_Or] :=
    (ZAnd[xp, #] & /@ orxp) // RemoveDuplicates
    
        (*TODO check that the structure is basically the same*)
ZAnd[xp_, andxp_And] :=
    With[ {fst = First[andxp], rst = Rest[andxp]},
    	(*Print[fst];*)
        Which[
            Head[fst] === Or,
            ReZAnd[xp,rst] /@ fst (*// RemoveDuplicates*),
            Head[fst]=== And,
            (*Print[andxp];*)
            ReZAnd[xp, rst, fst],
            True,
            ReZAnd[xp, rst, fst]
        ]
    ]
    
(*Operator form of ReZAnd*)    
ReZAnd[xp_,rst_] :=
    ReZAnd[xp,rst,#]&     
    
ReZAnd[xp_,rst_,fst_Equal] :=
    Module[ {fsol = First @ Solve @ fst // Quiet},
    	(*Print["pre right of ZAnd\n", rst];
    	Print["elements on the right hand side:", Length[rst]," sol: ", fsol];
    	Print["post right of ZAnd\n",ReplaceSolution[rst, fsol], "\n"];*)
        ZAnd[(xp /. fsol) && And @@ (fsol /. Rule -> Equal) // Simplify, ReplaceSolution[rst, fsol]]
    ]

(*TODO fix this if needed*)
ReZAnd[xp_, rst_, fst_And] :=
    (
    (*Print["got And: ", fst];*)
    ReZAnd[xp,rst] /@ fst 
    (*xp&&rst&&fst*)
    )  
ReZAnd[xp_,rst_,fst_Inequality] :=
    ZAnd[xp && fst, rst]
    
ReZAnd[xp_,rst_,fst_LessEqual] :=
    ZAnd[xp && fst, rst]

ReZAnd[xp_,rst_,fst_GreaterEqual] :=
    ZAnd[xp && fst, rst]
    
ReZAnd[xp_,rst_,fst_Less] :=
    ZAnd[xp && fst, rst]

ReZAnd[xp_,rst_,fst_Greater] :=
    ZAnd[xp && fst, rst]

ZAndEqual[rules_Association, andxp_And] :=
    With[ {fst = First[andxp], rst = Rest[andxp]},
    	Module[ {fsol = First @ Solve @ fst // Quiet},
    		(*Print[rules];*)
    		(*Print[andxp];
    		Print[fsol];*)
    		ZAndEqual[Join[rules/.fsol,Association[fsol]], ReplaceSolution[rst, fsol]]
    	]
    ]

ZAndEqual[rules_Association, xp_Equal] :=
	Module[ {fsol1 = First @ Solve @ xp // Quiet},
   		Join[rules/.fsol1, Association[fsol1]]
   	]
         
ReplaceSolution[True, sol_] :=
    True

ReplaceSolution[False, sol_] :=
    False

ReplaceSolution[rst_And, sol_] :=
    Module[{newrst},
    	newrst = rst /. sol;
    	If[Head[newrst] === And,
    		Reduce /@ (rst /. sol),
    		Reduce[rst /. sol]
    	]
    ]
    
    
ReplaceSolution[rst_, sol_] :=
    ((*Print["RS, rest: \n",rst,"\nrest with sol:\n",rst/.sol]; *)
    Simplify[rst /. sol])

RemoveDuplicates[xp_And] :=
    DeleteDuplicates[Sort[xp]];

RemoveDuplicates[xp_Or] :=
    DeleteDuplicates[Sort[xp]];
    
RemoveDuplicates[xp_] :=
    xp

NewReduce[s_Or] :=
    Reduce @ s

NewReduce[True] :=
    True

NewReduce[False] :=
    False
    
NewReduce[x_Equal] :=
(Print["here"];
    Reduce @ x)
    
NewReduce[s_Inequality] :=
    s
    
NewReduce[x_GreaterEqual] :=
    x

NewReduce[x_LessEqual] :=
    x

NewReduce[x_Inequality] :=
    x

NewReduce[system_And] :=
    Module[ {result,
    	groups = GroupBy[List @@ system, Head[#] === Or||Head[#]===Equal&],
    	subgroups},
    	subgroups = GroupBy[groups[True], Head[#]===Equal&];
    	result = ZAnd[And@@groups[False],And@@Lookup[subgroups,True,True]&&And@@subgroups[False]];
        (*result = ZAnd[Select[system, !((Head[#] === Or)||(Head[#]===Equal))&],Select[system, ((Head[#] === Or)||(Head[#]===Equal))&]];*)
        If[ result === False,
            False,
            result // DeleteDuplicates
        ]
    ]

NewCleanEqualities[{system_, rules_}] :=
	Module[{groups = GroupBy[List@@system, Head[#]===Equal&],
		newrules},
		newrules = ZAndEqual[rules, And@@groups[True]];
		{Simplify/@And@@(groups[False] /. newrules), newrules}
	]
	
NewCleanEqualities2[{system_, rules_}] :=
	Module[{groups = GroupBy[List@@system, Head[#]===Equal&],
		newrules},
		newrules = ZAndEqual2[And@@groups[False], Association@rules, And@@groups[True]];
		{Simplify/@And@@(groups[False] /. newrules), newrules}
	]

ZAndEqual2[xpo_, rules_Association, andxp_And] :=
    With[ {fst = First[andxp], rst = Rest[andxp]},
    	Module[ {fsol = First @ Solve @ fst // Quiet, newxpo, groups},
    		newxpo = ReplaceSolution[xpo, fsol];
    		groups = GroupBy[List @@ newxpo, Head[#] === Equal&];
    		(*Print[rules];*)
    		(*Print[andxp];*)
    		Print["head and ", fsol];
    		ZAndEqual2[And@@groups[False],Join[rules/.fsol,Association[fsol]], ReplaceSolution[rst, fsol]&&And @@ groups[False]]
    	]
    ]

ZAndEqual2[xpo_, rules_Association, xp_Equal] :=
	Module[ {fsol1 = First @ Solve @ xp // Quiet, newxpo, groups},
		newxpo = ReplaceSolution[xpo, fsol1];
   		groups = GroupBy[List @@ newxpo, Head[#] === Equal&];
   		Print["head equal ", groups[True]];
   		Join[rules/.fsol1, Association[fsol1]]
   	]
         
    
CleanEqualitiesOperator[Eqs_][system_] :=
    CleanEqualitiesOperator[Eqs][{system, {}}]

CleanEqualitiesOperator[Eqs_][{system_, rules_}] :=
    FixedPoint[EqEliminator[Eqs], {system,rules}]

EqEliminator[Eqs_][{system_, rules_List}] :=
    EqEliminator[Eqs][{system, Association[rules]}]

EqEliminator[Eqs_][{system_, rules_Association}] :=
    Module[ {EE, ON, newrules,groups},
    	Which[Head[system] === And,	
        (*separate equalities from the rest*)
        groups = GroupBy[List@@system, Head[#]===Equal&];
        EE = And@@Lookup[groups, True, Return[{system, rules},Module]];
        ON = And@@Lookup[groups, False, True],
        Print["And"];
        Head[system] === Equal,
        EE=system;
        ON = True,
        True,
        EE = True;
        ON = Simplify @ system;
    	];
        newrules = First @ Solve[EE] // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)
        newrules = Join[rules, Association @ newrules];
        {ON /. newrules, Simplify/@(newrules /. newrules)}
    ]

(*EqEliminator[Eqs_][{system_Equal, rules_Association}] :=
    Module[ {EE=system, newrules,ON},
            (*vars = Join[UOrder[Eqs],JOrder[Eqs],Values@Eqs["jtvars"]];(*Join[Eqs["uvars"]//Values,Eqs["jvars"]//Values,Eqs["jtvars"]//Values];*)*)
        ON =True;
        newrules = First @ Solve[EE]//Quiet;
        newrules = Join[rules, Association @ newrules];
        {ON/.newrules, Simplify/@(newrules/.newrules)}
    ]*)

(*EqEliminator[Eqs_][{True, rules_ }] :=
    {True, Simplify /@ rules}*)

(*EqEliminator[Eqs_][{system_Or, rules_ }] :=
    {Simplify @ system, Simplify /@ rules}*)
       
(*EqEliminator[Eqs_][{system_GreaterEqual, rules_ }] :=
    {Simplify @ system, Simplify /@ rules}*)
       
(*EqEliminator[Eqs_][{system_LessEqual, rules_ }] :=
    {Simplify @ system, Simplify /@ rules}*)
       
(*EqEliminator[Eqs_][{system_Inequality, rules_ }] :=
    {Simplify @ system, Simplify /@ rules}*)

(*EqEliminator[Eqs_][{False,rules_}] :=
    {False, rules}*)

FixedReduceX1[Eqs_Association][rules_] :=
    Module[ {nonlinear, auxsys, auxsol, error, TOL = 10^-10},
        nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], RoundValues[Eqs["Nrhs"]/. rules]}]);
        auxsol = First[Solve[nonlinear]//Quiet];
        auxsys = (Eqs["EqAllAll"] && nonlinear) /. auxsol;
        {auxsys, auxsol} = CleanEqualitiesOperator[Eqs][{auxsys, auxsol}];
        auxsys = NewReduce[auxsys];
        {auxsys, auxsol} = CleanEqualitiesOperator[Eqs][{auxsys, auxsol}];
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

CriticalCongestionSolver2[D2E_Association] :=
    Module[ {system, subsystem, rules, vars, sysclean, rulclean, originalsystem, originalrules,
    EqAllAll = Lookup[D2E, "EqAllAll", Print["No equations to solve."];
                                       Return[]], 
    EqCriticalCase = Lookup[D2E,"EqCriticalCase", Print["Critical case equations are missing."];
                                                  Return[]],
    InitRules = Lookup[D2E, "BoundaryRules", Print["Need boundary conditions"];],
    jts = Lookup[D2E, "jtvars"]//Values,
    js = JOrder[D2E],
    us = UOrder[D2E],
    nvars,
    newcrit,
    aux
},

    ]


CriticalCongestionSolver[D2E_Association] :=
    Module[ {system, subsystem, rules, vars, sysclean, rulclean, originalsystem, originalrules,
    EqAllAll = Lookup[D2E, "EqAllAll", Print["No equations to solve."];
                                       Return[]], 
    EqCriticalCase = Lookup[D2E,"EqCriticalCase", Print["Critical case equations are missing."];
                                                  Return[]],
    InitRules = Lookup[D2E, "BoundaryRules", Print["Need boundary conditions"];],
    jts = Lookup[D2E, "jtvars"]//Values,
    js = JOrder[D2E],
    us = UOrder[D2E],
    nvars,
    newcrit,
    aux
},
        (*Print["Number of variables: ", Length[Join[js, jts, us]]];*)
        system = EqCriticalCase && EqAllAll;
        rules = Association @ InitRules;
        
        {{newcrit,rules},aux} = FixedPoint[SimpleCrit,{{EqCriticalCase,rules}, Select[us, ! FreeQ[#][EqCriticalCase] &]}];
        originalsystem = system;
        originalrules = rules;
        (*system = Simplify /@ (system /. rules);*)
        (*system = (system /. rules);*)
        (*Print[system];*)
        system = Simplify/@(system/.rules);
        (*Print[system];*)
		{system, rules} = FixedPoint[SetUValuesStep[D2E], {system, rules}];
        {system, rules} = FixedPoint[SetJtValuesStep[D2E], {system, rules}];
        (*Print["before \"eliminate\" us\n",getEqual@system,"\n",Solve[getEqual@system]];*)
        {system, rules} = FixedPoint[SetJUValuesStep[D2E], {system, rules}];
        (*Print["before \"eliminate\" us\n",getEqual@system,"\n",Solve[getEqual@system]];*)
        {system, rules} = FixedPoint[SetJtUValuesStep[D2E], {system, rules}];
        (*Print["jts\n", getEqual@system];*)
        
        (*Print["before \"eliminate\" us\n",getEqual@system,"\n",Solve[getEqual@system]];*)
        Print[Length[rules]];
        	{{system, rules}, vars, nvars} = EliminateVarsSimplify[D2E][{{system, rules}, Join[us,js,jts], {}}];
        (*Print["\"eliminated\" us\n",getEqual@system,"\n",Solve[getEqual@system]];*)
        Print[Length[rules]];


        (*{{system, rules}, vars, nvars} = EliminateVars[D2E][{{system, rules}, js, {}}];*)
		(*Print["is u1 still here?\n", Simplify/@system];*)
        If[ AllTrue[FreeQ[#][system] & /@ Join[us, js] // DeleteDuplicates,TrueQ],
            Print["Done!"],
            vars = Select[Join[us, js], Function[var, !FreeQ[var][system] ] ];
            (*Print[Length[rules]];*)
            {{system, rules}, vars, nvars} = EliminateVarsSimplify[D2E][{{system, rules}, vars, {}}];
(*            {system, rules} = CleanEqualitiesOperator[D2E][{Simplify@system, rules}]*)
        ];
        {RemoveDuplicates[system/.rules], rules} /. rules
    ]

CriticalCongestionStep[{sys_Or, rul_}] :=
    {BooleanConvert[sys,"CNF"], rul}

CriticalCongestionStep[{False, rul_}] :=
    {False, rul}

CriticalCongestionStep[{True, rul_}] :=
    {True, rul}

CriticalCongestionStep[{sys_,rul_}] :=
    (Print["CCS: Removed variables: ",Length[rul]];
     CleanEqualitiesOperator[D2E][{NewReduce[BooleanConvert[sys,"CNF"](*sys*)], rul}])

EliminateVarsStep[D2E_][{{system_, rules_}, {}, {}}] :=
    {{system, rules}, {}, {}}

EliminateVarsSimplifyStep[D2E_][{{system_, rules_}, {}, {}}] :=
    {{system, rules}, {}, {}}
 
 
 (*Coloque o sistema original com as regras substituidas*)
EliminateVarsStep[D2E_][{{system_, rules_}, us_, persistus_}] :=
    Module[ {var, subsys, position , 
    	newsys = system, newrules = Association @ rules, 
    	newus, subsyscomplement, newpersistus = persistus,
    	groups},
        var = SelectFirst[us, ! FreeQ[#][newsys] &];
        If[ Head[var] === Missing,
            Return[{{newsys, newrules}, {}, persistus}]
        ];
        (*Print["Selecting expressions with ", var];*)
        position = First @ FirstPosition[us, var];
        newus = Drop[us, position];
        groups = GroupBy[List@@newsys, FreeQ[var]];
        (*Print["EliminateVarsStep: \n"];*)
        subsys = Simplify/@And@@groups[False];
        (*Print["before cleaneq:\n",subsys,"\ngetequal\n",getEqual@subsys, "\n",newrules];*)
        {subsys, newrules} = CleanEqualitiesOperator[D2E][{subsys, newrules}];
        (*Print["after cleaneq:\n",subsys, "\n",newrules];*)
        newsys = newsys /. newrules;
        If[ subsys =!= True,
            newsys = subsys && And@@groups[True];
            If[ !FreeQ[var][subsys],
                    (*Print[var];*)
                newpersistus = Join[newpersistus,{var}]
            (*Print["Removed ", var];*)
                ];
        ];
        {{Simplify/@(newsys/.newrules), Simplify/@(newrules/.newrules)}, newus, newpersistus}
    ]
    
EliminateVarsSimplifyStep[D2E_][{{system_, rules_}, us_, persistus_}] :=
    Module[ {var, subsys, position , newsys = system, newrules = Association @ rules, newus, subsyscomplement, newpersistus = persistus},
        var = SelectFirst[us, ! FreeQ[#][newsys] &];
        If[ Head[var] === Missing,
            Return[{{newsys, newrules}, {}, persistus}]
        ];
        (*Print["Selecting expressions with ", var];*)
        position = First @ FirstPosition[us, var];
        newus = Drop[us, position];
        subsys = Select[newsys, Function[x, !FreeQ[var][x]]];
        (*Print[subsys];*)
        subsys = Simplify @ subsys;
        {subsys, newrules} = CleanEqualitiesOperator[D2E][{subsys, newrules}];
        (*Print[subsys];*)
        newrules = Simplify /@ newrules;
        newsys = newsys /. newrules;
        If[ subsys =!= True,
            subsyscomplement = Select[newsys, Function[x, FreeQ[var][x]]];
            newsys = subsys && subsyscomplement;
            If[ !FreeQ[var][subsys],
                    (*Print[var];*)
                newpersistus = Join[newpersistus,{var}]
            (*Print["Removed ", var];*)
                ];
        ];
        {{newsys, newrules}, newus, newpersistus}
    ]
EliminateVarsSimplify[D2E_][{{system_, rules_}, us_, persistus_}] :=
    FixedPoint[EliminateVarsSimplifyStep[D2E], {{system, rules}, us, persistus}]     
EliminateVars[D2E_][{{system_, rules_}, us_, persistus_}] :=
    FixedPoint[EliminateVarsStep[D2E], {{system, rules}, us, persistus}](*(Length[#1[[1,2]]] === Length[#2[[1,2]]]&)]*)

SimpleCrit[{{True, rules_}, critus_}] := {{True, rules}, critus}
SimpleCrit[{{syscrit_, rules_}, {}}] :=
 Module[{subsol, newcrit = Simplify /@ syscrit},
  subsol = First@Solve[newcrit];
  {{newcrit /. subsol, Join[rules, Association[subsol]] /. subsol}, {}}
  ]

SimpleCrit[{{syscrit_, rules_}, critus_}] :=
 Module[{var, subsys, subsol, newcrit = Simplify /@ syscrit},
  var = First[critus];
  subsys = Select[newcrit, Function[xp, ! FreeQ[var][xp]]];
  If[subsys === True,
   Return[{{newcrit, rules}, Rest[critus]}]
   ];
  subsol = First@Solve[subsys[[1]], var];
  {{syscrit /. subsol, Join[rules, Association[subsol]] /. subsol}, 
   Rest[critus]}
  ]
 
End[]