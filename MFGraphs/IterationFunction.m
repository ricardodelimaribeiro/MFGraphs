(* ::Package:: *)

(* Wolfram Language package *)

CriticalCongestionSolver2::usage = 
"CriticalCongestionSolver[eq_Association] returns the critical congestion solution."

CleanEqualitiesOperator::usage = 
"CleanEqualitiesOperator[Eqs][{system,rules}] returns an equivalent pair of system and rules without \"loose\" equalities. 
The equalities are solved and substituted into the system and added to rules."

CleanEqualities::usage = 
"CleanEqualitiesOperator[Eqs][{system,rules}] returns an equivalent pair of system and rules without \"loose\" equalities. 
The equalities are solved and substituted into the system and added to rules."

EqEliminator::usage = 
"EqEliminator[Eqs][{system, rules}] extracts the equalities of the system, solves and appends them to rules, and returns a new pair {system, rules} /. rules. ";

NewReduce::usage = 
"NewReduce[system] reduces the system using ZAnd"

ZAnd::usage =
"ZAnd[system, orsys] is a recursive function to reduce system&&orsys by checking the feasibility of the alternatives in orsys."

EliminateVars::usage = 
"EliminateVars[Eqs][{{system, rules}, us, persistus}]"

EliminateVarsSimplify::usage = 
""

Begin["`Private`"]
EliminateVarsStep::usage = 
""

RemoveDuplicates::usage =
"RemoveDuplicates[exp] Sort and DeleteDuplicates"

CriticalCongestionStep::usage =
""

ReZAnd::usage =
"ReZAnd[] "

Solver::usage =
"Solver[Eqs_Association][alpha] solves the network problem with the given alpha."

FixedReduceX1::usage = 
"FixedReduceX1[MFGEquations][rules] gives the solution for the system when the RHS of the nonlinear equations have rules substituted in them.
	MFGEquations should have the Keys EqAllAll, BoundaryRules, Nlhs, Nrhs, TOL";(*TODO review the keys when critical solver is done...*)

CleanAndReplace::usage =
"CleanAndReplace[{system,rules}] returns True if the system is approximately (difference in each equation is less than 10^(-10) ) solved by rules.";

CleanAndReplace[{system_,rules_}] :=
    (system /. Equal -> (zob[#1-#2]&)/. rules // Chop(*[#,10^(-12)]&*) ) /. zob -> (# == 0.&)
 
ZAnd[_, False] :=
    (
     False)

ZAnd[xp_, leq_LessEqual] :=
    Simplify[xp && leq]

ZAnd[xp_, geq_GreaterEqual] :=
    Simplify[xp && geq]

ZAnd[xp_, ineq_Inequality] :=
    Simplify[xp && ineq]

ZAnd[False, _] :=
    (
     False)

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
            Head[fst]=== And,
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
        ZAnd[(xp /. fsol) && And @@ (fsol /. Rule -> Equal) // Simplify, ReplaceSolution[rst, fsol]]
    ]

ReZAnd[xp_, rst_, fst_And] :=
    ReZAnd[xp,rst] /@ fst
    
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

ReplaceSolution[True, sol_] :=
    True

ReplaceSolution[False, sol_] :=
    False

ReplaceSolution[rst_And, sol_] :=
    Module[{newrst},
    	newrst = rst /. sol;
    	If[Head[newrst] === And,
    		Reduce[#,Reals]& /@ (rst /. sol),
    		Reduce[rst /. sol,Reals]
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
    Reduce[s,Reals]

NewReduce[True] :=
    True

NewReduce[False] :=
    False
    
NewReduce[x_Equal] :=
(Print["here"];
    Reduce[x,Reals])
    
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
    	subgroups,
    	sorted},
    	subgroups = GroupBy[groups[True], Head[#]===Equal&];
    	If[Head[subgroups[False]]===Missing,
    		result = ZAnd[And@@groups[False],And@@Lookup[subgroups,True,True]],
    		sorted = SortBy[And@@Lookup[subgroups,False, True],Simplify`SimplifyCount];
    		result = ZAnd[And@@groups[False],(And@@Lookup[subgroups,True,True])&&sorted]
    	];
        (*result = ZAnd[Select[system, !((Head[#] === Or)||(Head[#]===Equal))&],Select[system, ((Head[#] === Or)||(Head[#]===Equal))&]];*)
        If[ result =!= False,
            result = result // DeleteDuplicates
        ];
        result
    ]

CleanEqualitiesOperator[Eqs_][system_] :=
    CleanEqualitiesOperator[Eqs][{system, {}}]

CleanEqualitiesOperator[Eqs_][{system_, rules_}] :=
    FixedPoint[EqEliminator[Eqs], {system /. rules, rules}]

CleanEqualities[system_]:=
	CleanEqualities[{system,{}}]
CleanEqualities[{system_, rules_}] :=
	CleanEqualitiesOperator[<||>][{system, rules}]

CleanEqualities[{system_List, rules_}] :=
	Module[{systems, ruless, aux},
	aux = CleanEqualitiesOperator[<||>][{#, rules}] &/@ system;
	(*aux =First @ aux;*)
	systems = First /@ aux;
	systems = And @@ systems;
	ruless = Last /@ aux;
	ruless = Join @@ ruless;
	systems = systems /. ruless;
	ruless = ruless /. ruless;
	{systems, ruless}
	]

EqEliminator[Eqs_][{system_, rules_List}] :=
    EqEliminator[Eqs][{system, Association[rules]}]

EqEliminator[Eqs_][{system_, rules_Association}] :=
    Module[ {EE, ON, newrules,groups},
    	Which[Head[system] === And,	
        (*separate equalities from the rest*)
        groups = GroupBy[List@@system, Head[#]===Equal&];
        EE = And@@Lookup[groups, True, Return[{system, rules},Module]];
        ON = And@@Lookup[groups, False, True],
        Head[system] === Equal,
        EE = system;
        ON = True,
        True,
        EE = True;
        ON = Simplify @ system;
    	];
        newrules = First @ Solve[EE] // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)
        newrules = Join[rules, Association @ newrules];
        {ON /. newrules, Expand/@(newrules /. newrules)}
    ]

EqEliminator[{system_, rules_}]:= EqEliminator[<||>][{system, rules}]

CriticalCongestionSolver2[Eqs_Association] :=
    Module[ {system, 
    AllIneqs,
	AllOrs,
    rules,
    js = Values @ Lookup[Eqs, "jvars"],
    aux,
    rvars,
    uvalues,
    jsys,
    jrules,
    us = Values[Eqs["uvars"]],
    newus,
    newjs
    },
    AllIneqs = Lookup[Eqs, "AllIneq", Print["No inequalities to solve."];
                                       Return["a"]];
	AllOrs = Lookup[Eqs, "AllOr", Print["No alternatives to solve. \n",AllOrs];
                                       Return["b"]];
    rules = Lookup[Eqs,"RulesCriticalCase", Print["Critical case rules are missing."];
                                                  Return[]];
    
    (*Replace critical congestion rules!*)
    AllOrs = AllOrs /. rules;
    AllIneqs = AllIneqs /. rules;
    (*try to simplify bit by bit:*)
    system = AllIneqs&&AllOrs;
    rvars = Variables[Join[us,js]/.rules];
    newus = Intersection[us,rvars];
    Print["EliminateVarsSimplify for the us"];
    {{system, rules}, aux, newus} = EliminateVarsSimplify[Eqs][{{system, rules}, newus, {}}];
    rules = Expand/@rules;
    If[Variables[us/.rules] === {},
    	Print["Finished with the us!"];
    	If[Variables[js/.rules] === {},
    		Print["Done!"],
    		Print["Finish with the js"];
    		uvalues = AssociationThread[us, us /. rules];
    		jsys = (Eqs["EqCriticalCase"] /. uvalues) && Eqs["EqPosJs"] && Eqs["EqCurrentCompCon"];
    		(**Retrieve js values already defined*)
    		jrules = Select[AssociationThread[js,js/.rules], NumericQ];
    		{jsys, jrules} = CleanEqualities[{jsys,jrules}];
    		{{jsys, jrules}, aux, newjs} = EliminateVarsSimplify[Eqs][{{jsys, jrules}, js}];
    		{jsys, jrules} = CleanEqualities[{jsys,jrules}];
    		
  			AssociateTo[uvalues, jrules];
    	],  
    	Print["Do what?"];
    	
    ];
	uvalues
]

EliminateVarsStep[Eqs_][{{system_, rules_}, {}, {}}] :=
    {{system, rules}, {}, {}}
 
(*Coloque o sistema original com as regras substituidas*)
EliminateVarsStep[Eqs_][{{system_, rules_}, us_, persistus_}] :=
    Module[ {var, subsys, position , 
    	newsys = system, newrules = Association @ rules, 
    	newus, newpersistus = persistus,
    	groups},
        var = SelectFirst[us, ! FreeQ[newsys, #] &];
        If[ Head[var] === Missing,
            Return[{{newsys, newrules}, {}, persistus}]
        ];
        position = First @ FirstPosition[us, var];
        newus = Drop[us, position];
        groups = GroupBy[List@@newsys, FreeQ[#, var]&];
        subsys = Simplify@And@@groups[False];
        {subsys, newrules} = CleanEqualitiesOperator[Eqs][{subsys, newrules}];
        newsys = newsys /. newrules;
        If[ subsys =!= True,
            newsys = subsys && And@@groups[True];
            If[ !FreeQ[subsys, var],
                newpersistus = Join[newpersistus,{var}]
                ];
        ];
        {{Simplify/@(newsys/.newrules), Simplify/@(newrules/.newrules)}, newus, newpersistus}
    ]
    
EliminateVarsSimplifyStep[Eqs_][{{system_, rules_}, {}, {}}] :=
    {{system, rules}, {}, {}}
 
EliminateVarsSimplifyStep[Eqs_][{{True, rules_}, us_, persistus_}] :=
    {{True, rules}, us, persistus}
 
EliminateVarsSimplifyStep[Eqs_][{{system_, rules_}, us_, persistus_}] :=
    Module[ 
    	{var, subsys, position , newsys = system, 
    	newrules = Association @ rules, newus, 
    	newpersistus = persistus, groups
    	},
    	groups = GroupBy[List @@ system, Head[#] === Or &];
        var = SelectFirst[us, ! FreeQ[And @@ Lookup[groups, True,True], #] &];
        If[ Head[var] === Missing,
            Return[{{newsys, newrules}, {}, persistus}]
        ];
        position = First @ FirstPosition[us, var];
        newus = Drop[us, position];
        If[newus == {},
        	subsys = (*And @@ Lookup[groups, False, True] &&*) Select[newsys, Function[x, !FreeQ[x, var]]]
        	,
        	subsys = (*And @@ Lookup[groups, False, True] &&*) Select[newsys, Function[x, !(FreeQ[x, var] && FreeQ[x, First @ newus])]]
        ];
        (*Print["EliminateVarsSimplifyStep: Reducing ..."];*)
        subsys = Reduce[subsys];
        subsys = Reduce[subsys, Reals];
        (*Print["EliminateVarsSimplifyStep: Clean equalities for system of ", var];*)
        {subsys, newrules} = CleanEqualitiesOperator[Eqs][{subsys, newrules}];
        newsys = newsys /. newrules;
        {{newsys, newrules}, newus, newpersistus}
    ]
    
EliminateVarsSimplify[Eqs_][{{system_, rules_}, us_}] :=
    FixedPoint[EliminateVarsSimplifyStep[Eqs], {{system /. rules, rules}, us, {}}]     

EliminateVarsSimplify[Eqs_][{{system_, rules_}, us_, persistus_}] :=
    FixedPoint[EliminateVarsSimplifyStep[Eqs], {{system /. rules, rules}, us, persistus}]     

EliminateVars[Eqs_][{{system_, rules_}, us_, persistus_}] :=
    FixedPoint[EliminateVarsStep[Eqs], {{system /. rules, rules}, us, persistus}](*(Length[#1[[1,2]]] === Length[#2[[1,2]]]&)]*)

EliminateVars[Eqs_][{{system_, rules_}, us_}] := 
	EliminateVars[Eqs][{{system, rules}, us, {}}]
 
End[]
(**********************************************************************)
(*
CriticalCongestionSolver::usage = 
"CriticalCongestionSolver[eq_Association] returns the critical congestion solution."


CriticalCongestionSolver[Eqs_Association] :=
    Module[ {system, rules, vars, originalsystem, originalrules,
    EqAllAll = Lookup[Eqs, "EqAllAll", Print["No equations to solve."];
                                       Return[]], 
    EqCriticalCase = Lookup[Eqs,"EqCriticalCase", Print["Critical case equations are missing."];
                                                  Return[]],
    InitRules = Lookup[Eqs, "BoundaryRules", Print["Need boundary conditions"];],
    jts = Lookup[Eqs, "jtvars"]//Values,
    js = JOrder[Eqs],
    us = UOrder[Eqs],
    nvars,
    newcrit,
    aux
},
        (*Print["Number of variables: ", Length[Join[js, jts, us]]];*)
        system = EqCriticalCase && EqAllAll;
        rules = Association @ InitRules;
        
        {{newcrit,rules},aux} = FixedPoint[SimpleCrit,{{EqCriticalCase,rules}, Select[us, ! FreeQ[EqCriticalCase, #] &]}];
        originalsystem = system;
        originalrules = rules;
        (*system = Simplify /@ (system /. rules);*)
        (*system = (system /. rules);*)
        (*Print[system];*)
        system = Simplify/@(system/.rules);
        (*Print[system];*)
		{system, rules} = FixedPoint[SetUValuesStep[Eqs], {system, rules}];
        {system, rules} = FixedPoint[SetJtValuesStep[Eqs], {system, rules}];
        (*Print["before \"eliminate\" us\n",getEqual@system,"\n",Solve[getEqual@system]];*)
        {system, rules} = FixedPoint[SetJUValuesStep[Eqs], {system, rules}];
        (*Print["before \"eliminate\" us\n",getEqual@system,"\n",Solve[getEqual@system]];*)
        {system, rules} = FixedPoint[SetJtUValuesStep[Eqs], {system, rules}];
        (*Print["jts\n", getEqual@system];*)
        
        (*Print["before \"eliminate\" us\n",getEqual@system,"\n",Solve[getEqual@system]];*)
        Print[Length[rules]];
        	{{system, rules}, vars, nvars} = EliminateVarsSimplify[Eqs][{{system, rules}, Join[us,js,jts], {}}];
        (*Print["\"eliminated\" us\n",getEqual@system,"\n",Solve[getEqual@system]];*)
        Print[Length[rules]];


        (*{{system, rules}, vars, nvars} = EliminateVars[Eqs][{{system, rules}, js, {}}];*)
		(*Print["is u1 still here?\n", Simplify/@system];*)
        If[ AllTrue[FreeQ[system, #] & /@ Join[us, js] // DeleteDuplicates,TrueQ],
            Print["Done!"],
            vars = Select[Join[us, js], Function[var, !FreeQ[system, var] ] ];
            (*Print[Length[rules]];*)
            {{system, rules}, vars, nvars} = EliminateVarsSimplify[Eqs][{{system, rules}, vars, {}}];
(*            {system, rules} = CleanEqualitiesOperator[Eqs][{Simplify@system, rules}]*)
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
     CleanEqualitiesOperator[Eqs][{NewReduce[BooleanConvert[sys,"CNF"](*sys*)], rul}])




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


*)
(**********************************************************************)



(*
CriticalBundle[Data_Association] :=
    Module[ {d2e},
        d2e = D2E[Data];
        {d2e,CriticalCongestionSolver[d2e]}
    ]


ZAndEqual::usage =
"ZAndEqual[rules, xp] returns an association with the solution of the equalities"

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
   	
NewCleanEqualities::usage =(*this strategy is too slow! maybe we can repair...*)
"NewCleanEqualities[{system, rules}] returns the non-equalities and the solution to the complementary linear system."

SimpleCrit::usage=
"SimpleCrit[{{syscrit, rules}, critus}] solves the critical congestion equations for the us."

SimpleCrit[{{True, rules_}, critus_}] := {{True, rules}, critus}
SimpleCrit[{{syscrit_, rules_}, {}}] :=
 Module[{subsol, newcrit = Simplify /@ syscrit},
  subsol = First@Solve[newcrit];
  {{newcrit /. subsol, Join[rules, Association[subsol]] /. subsol}, {}}
  ]

SimpleCrit[{{syscrit_, rules_}, critus_}] :=
 Module[{var, subsys, subsol, newcrit = Simplify /@ syscrit},
  var = First[critus];
  subsys = Select[newcrit, Function[xp, ! FreeQ[xp, var]]];
  If[subsys === True,
   Return[{{newcrit, rules}, Rest[critus]}]
   ];
  subsol = First@Solve[subsys[[1]], var];
  {{syscrit /. subsol, Join[rules, Association[subsol]] /. subsol}, 
   Rest[critus]}
  ]

*)         
