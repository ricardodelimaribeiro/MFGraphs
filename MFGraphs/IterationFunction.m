(* ::Package:: *)

(* Wolfram Language package *)

CriticalCongestionSolver2::usage = 
"CriticalCongestionSolver2[eq_Association] returns the critical congestion solution."

CriticalCongestionSolver::usage =
"CriticalCongestionSolver[eq_Association] returns the critical congestion solution. Developing solver for non-zero switching costs."

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

EliminateVarsSimplify::usage = 
"EliminateVarsSimplify[Eqs_][{{system, rules}, us, persistus}] reduces the subsystems with us"

EliminateVars::usage = 
"EliminateVars[Eqs_][{{system, rules}, us, persistus}] reduces the subsystems with us"

FixedReduce2::usage =
"FixedReduce2[Eqs_Association][rules_] is the step for the non-linear solver"

FixedReduce3::usage =
"FixedReduce3[Eqs_Association][rules_] is the step for the non-linear solver.
This one defines the u's from the approximating currents, then redefines the currents."

Begin["`Private`"]
RemoveDuplicates::usage =
"RemoveDuplicates[exp] Sort and DeleteDuplicates"

ReZAnd::usage =
"ReZAnd[] "

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
    Module[ {newrst},
        newrst = rst /. sol;
        If[ Head[newrst] === And,
            Reduce[#,Reals]& /@ (rst /. sol),
            Reduce[rst /. sol,Reals]
        ]
    ]
    
ReplaceSolution[rst_, sol_] :=
    Simplify[rst /. sol]

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
    
(*TODO NewReduce: is there a problem here?*)
NewReduce[system_And] :=
    Module[ {result,
        groups = GroupBy[List @@ system, Head[#] === Or||Head[#]===Equal&],
        subgroups,
        sorted},
        subgroups = GroupBy[groups[True], Head[#]===Equal&];
        If[ Head[subgroups[False]]===Missing,
            result = ZAnd[And@@Lookup[groups,False, True],And@@Lookup[subgroups,True,True]],
            sorted = SortBy[And@@Lookup[subgroups,False, True],Simplify`SimplifyCount];
            result = ZAnd[And@@Lookup[groups, False, True],(And@@Lookup[subgroups,True,True])&&sorted]
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

CleanEqualities[system_] :=
    CleanEqualities[{system,{}}]
    
CleanEqualities[{system_, rules_}] :=
    CleanEqualitiesOperator[<||>][{system, rules}]

CleanEqualities[{system_List, rules_}] :=
    Module[ {systems, ruless, aux},
        aux = CleanEqualitiesOperator[<||>][{#, rules}] &/@ system;
        systems = First /@ aux;
        systems = And @@ systems;
        ruless = Last /@ aux;
        ruless = Join @@ ruless;
        systems = systems /. ruless;
        ruless = ruless /. ruless;
        {systems, ruless}
    ]

EqEliminator[{system_, rules_}] :=
    EqEliminator[<||>][{system, rules}]

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
    
    
FixedReduce2[Eqs_Association][<||>] :=
    CriticalCongestionSolver2[Eqs]   

FixedReduce2[Eqs_Association][{}] :=
    CriticalCongestionSolver2[Eqs]   
    
FixedReduce2[Eqs_Association][approxrules_Association] :=
    Module[ {(*auxsys, auxsol, error, TOL = 10^-10,*) system, 
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
        newus, RuleNonCritical,
        newjs
        },
        (*nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], RoundValues[Eqs["Nrhs"]/. approxrules]}]);
        nonlinear = nonlinear /. Eqs["InitRules"];*)
        (*Lets use "RuleNonCritical" instead.*)
(*        Print[Eqs["Nrhs"]];
        Print["nonlinear ewqs\n",nonlinear];
        Print["replacement rules: \n", Eqs["costpluscurrents"]/.approxrules];*)
        AllIneqs = Lookup[Eqs, "AllIneq", Print["No inequalities to solve."];
                                          Return["a"]];
        AllOrs = Lookup[Eqs, "AllOr", Print["No alternatives to solve. \n",AllOrs];
                                      Return["b"]];
        RuleNonCritical = Eqs["RuleNonCritical"]/.Eqs["InitRules"];
        rules = Expand /@ (RuleNonCritical/.RoundValues[Eqs["costpluscurrents"]/.approxrules]);
        Print["non critical rules with j values:\n",rules];
        (*Replace critical congestion rules!*)
        AllOrs = AllOrs /. rules;
        AllIneqs = AllIneqs /. rules;
        
        (*try to simplify bit by bit:*)
        system = AllIneqs&&AllOrs;
        rvars = Variables[Join[us,js]/.rules];
        Print["FR2: EliminateVarsSimplify for the us"];
        {{system, rules}, aux, newus} = EliminateVarsSimplify[Eqs][{{system, rules}, (*new*)us, {}}];
        rules = Expand/@rules;
        If[ Variables[us/.rules] === {},
            Print["FR2: Finished with the us!"],
            Print["FR2: Not finished with the us!"];
            {{system, rules}, aux, newus} = EliminateVars[Eqs][{{system, rules}, (*new*)us, {}}];
            Print["FR2: Reducing further (NewReduce)..."(*, system*)];
            system = NewReduce[system];
            Print["FR2: Simplifying..."];
            system = Simplify @ system;
            {system, rules} = CleanEqualities[{system, rules}];
            uvalues = AssociationThread[us, us /. rules];
            jrules = AssociationThread[js, js /. rules];
            AssociateTo[uvalues, jrules];
            If[ Variables[us/.uvalues] =!= {},
                Print["FR2: (Maybe the) system is not feasible!\n
            \tCheck the paper: The current method for stationary mean-field games on networks"]
            ];
        ];
        If[ Variables[js/.rules] === {},
            Print["FR2: Done!\n"];
            Return[AssociationThread[Join[us, js], Join[us, js] /. rules]],
            Print["FR2: Finish with the js"];
            uvalues = AssociationThread[us, us /. rules];
            (*TODO Fix this! Use rules or correct the equations! not the criticalcase!*)
            jsys = (*(Eqs["EqCriticalCase"] /. uvalues) &&*) Eqs["EqPosJs"] && Eqs["EqCurrentCompCon"];
            (**Retrieve js values already defined*)
            (*Keeping the order for the js Association*)
            jrules = AssociationThread[js, js/.Select[AssociationThread[js, js/.rules], NumericQ]];
            Print["jrules: \n", jrules];
            {jsys, jrules} = CleanEqualities[{jsys, jrules}];
            {{jsys, jrules}, aux, newjs} = EliminateVarsSimplify[Eqs][{{jsys, jrules}, js}];
            {jsys, jrules} = CleanEqualities[{jsys,jrules}];
            AssociateTo[uvalues, jrules];
        ];
        Print["returned uvalues: \n",Values@uvalues//N];
        uvalues
    ]
 (*    Pause[30];       

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
        (*auxsol*)
 *)
 
 
 
CriticalCongestionSolver2[Eqs_Association] :=
    Module[ {system, AllIneqs,AllOrs,
    rules, js = Values @ Lookup[Eqs, "jvars"],
    aux, uvalues,
    jsys, jrules, us = Values[Eqs["uvars"]],
    newus, newjs, orineqs
    },
        AllIneqs = Lookup[Eqs, "AllIneq", Print["CCS2: No inequalities to solve."];
                                          Return["a"]];
        AllOrs = Lookup[Eqs, "AllOr", Print["CCS2: No alternatives to solve. \n",AllOrs];
                                      Return["CCS2: b"]];
        rules = Lookup[Eqs,"RulesCriticalCase", Print["CCS2: Critical case rules are missing."];
                                                Return[]];
        (*Replace critical congestion rules!*)
        AllOrs = AllOrs /. rules;
        AllIneqs = Simplify[AllIneqs /. rules];
        (*AllIneqs = AllIneqs /. rules;*)
        (*Print[Simplify[#, AllIneqs]& /@ AllOrs];*)
        (*AllOrs = Reduce[#, Reals]& /@ AllOrs;*)(*should remove non-equalities: !=*)
        Print[AllIneqs];
        AllOrs = Simplify[#, AllIneqs]& /@ AllOrs;
        Print[AllOrs];
        {system, rules} = CleanEqualities[{AllIneqs && AllOrs, rules}];
        Print[system];
        (*orineqs = GroupBy[List @@ AllOrs, Head[#]===Or&];
        Print[orineqs];
        AllOrs = And @@ Lookup[orineqs, True, True];
        AllIneqs = AllIneqs && And @@ Lookup[orineqs, False, True];
        AllIneqs = Simplify[AllIneqs];
        AllOrs = Simplify[#, AllIneqs]& /@ AllOrs;
        Print[AllOrs];
        AllOrs = BooleanConvert[AllOrs, "CNF"];
        (*Print["2\n",AllIneqs];
        Print["3\n",AllOrs];*)
        system = AllIneqs&&AllOrs;*)
        (*system = FullSimplify/@system;*)
        
        Print["CCS2: EliminateVarsSimplify for the us"];
        {{system, rules}, aux, newus} = EliminateVarsSimplify[Eqs][{{system, rules}, us, {}}];
        rules = Expand /@ rules;
        (*Print[system];
        Print[rules];*)
        If[ Variables[us /. rules] === {},
            Print["CCS2: Finished with the us!"],
            Print["CCS2: Not finished with the us!"];
            system = BooleanConvert[system, "CNF"];
            (*Print[system];*)
            system = BooleanConvert[#, "CNF"]& /@ system;
	        (*Print[system];*)
            {{system, rules}, aux, newus} = EliminateVars[Eqs][{{system, rules}, us, {}}];
            (*Print[system];*)
            system = Simplify /@ system;
            Print["CCS2: Reducing further (NewReduce)..."];
            system = NewReduce[system];
            Print["CCS2: Simplifying..."];
            system = Simplify @ system;
            {system, rules} = CleanEqualities[{system, rules}];
            uvalues = AssociationThread[us, us /. rules];
            jrules = AssociationThread[js, js /. rules];
            AssociateTo[uvalues, jrules];
            If[ Variables[us/.uvalues] =!= {},
                Print["CCS2: (Maybe the) system is not feasible!\n
            \tCheck the paper: The current method for stationary mean-field games on networks"],
                Print["CCS2: Finished with the us!"];
            ];
        ];
        If[ Variables[js/.rules] === {},
            Print["CCS2: Finished with the js!"];
            Return[AssociationThread[Join[us, js], Join[us, js] /. rules]],
            Print["CCS2: Now, finish with the js"];
            uvalues = AssociationThread[us, us /. rules];
            jsys = (Eqs["EqCriticalCase"] /. uvalues) && Eqs["EqPosJs"] && Eqs["EqCurrentCompCon"];
            (**Retrieve js values already defined*)
            (*Keeping the order for the js Association*)
            jrules = AssociationThread[js, js/.Select[AssociationThread[js, js/.rules], NumericQ]];
            {jsys, jrules} = CleanEqualities[{jsys, jrules}];
            {{jsys, jrules}, aux, newjs} = EliminateVarsSimplify[Eqs][{{jsys, jrules}, js}];
            {jsys, jrules} = CleanEqualities[{jsys,jrules}];
            AssociateTo[uvalues, jrules];
        ];
        uvalues
    ]

EliminateVarsSimplifyStep[Eqs_][{{system_, rules_}, us_, persistus_}] :=
    Module[ {subsys, newsys, 
        newrules, newus, allus
        },
        allus = Values[Eqs["uvars"]];
        newus = Drop[us, UpTo[1]];
        subsys = Select[system, Function[x, !(FreeQ[x, #])]]&/@us;
        subsys = RemoveDuplicates[DeleteCases[subsys,True]];
        (*Print["subsys list no duplicates\n",subsys];*)
        subsys = Reduce /@ subsys;(*Imaginary numbers may artificially appear: I * Im[jt556] ...*)
        subsys = Simplify[#, Eqs["EqPosJts"] && Eqs["EqPosJs"]]& /@ subsys;
        (*subsys = BooleanConvert[#,"CNF"]& /@ subsys;*)
        subsys = And @@ subsys;
        {subsys, newrules} = CleanEqualitiesOperator[Eqs][{subsys, rules}];
        newsys = system /. newrules;
        newsys = Reduce[#, Reals]& /@ newsys;
        (*newsys = FullSimplify /@ newsys;*)
        {{newsys, newrules}, us, Intersection[allus, getVar[newsys]]}
    ]
    
EliminateVarsSimplifyStep[Eqs_][{{system_, rules_}, {}, {}}] :=
    {{system, rules}, {}, {}}
 
EliminateVarsSimplifyStep[Eqs_][{{True, rules_}, us_, persistus_}] :=
    {{True, rules}, us, persistus}
 
EliminateVarsSimplify[Eqs_][{{system_, rules_}, us_}] :=
    EliminateVarsSimplify[Eqs][{{system /. rules, rules}, us, {}}]     

EliminateVarsSimplify[Eqs_][{{system_, rules_}, us_, persistus_}] :=
    FixedPoint[EliminateVarsSimplifyStep[Eqs], {{system /. rules, rules}, us, persistus}]     

CriticalCongestionSolver[Eqs_Association] :=
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
        {{system, rules}, aux, newus} = EliminateVars[Eqs][{{system, rules}, newus, {}}];
        rules = Expand/@rules;
        If[ Variables[us/.rules] === {},
            Print["Finished with the us!"],
            Print["(Maybe the) system is not feasible!\n
            \tCheck the paper: The current method for stationary mean-field games on networks"];
            uvalues = AssociationThread[us, us /. rules];
            Print[uvalues];
            uvalues = Select[uvalues, NumericQ];
            Print[uvalues];
            jsys = (Eqs["EqCriticalCase"] /. uvalues) && Eqs["EqPosJs"] && Eqs["EqCurrentCompCon"];
            jrules = AssociationThread[js, js/.Select[AssociationThread[js, js/.rules], NumericQ]];
            jrules = Join[uvalues, jrules];
            Print[jrules];
            {jsys, jrules} = CleanEqualities[{jsys, jrules}];
            Print[rules];
            Print[Expand/@(rules/.(jrules/.rules))];
            Print[jsys];
        ];
        If[ Variables[js/.rules] === {},
            Print["Done!"];
            Return[AssociationThread[Join[us, js], Join[us, js] /. rules]],
            Print["Finish with the js"];
            uvalues = AssociationThread[us, us /. rules];
            Print[uvalues];
            jsys = (Eqs["EqCriticalCase"] /. uvalues) && Eqs["EqPosJs"] && Eqs["EqCurrentCompCon"];
            (**Retrieve js values already defined*)
            (*Keeping the order for the js Association*)
            jrules = AssociationThread[js, js/.Select[AssociationThread[js, js/.rules], NumericQ]];
            {jsys, jrules} = CleanEqualities[{jsys, jrules}];
            {{jsys, jrules}, aux, newjs} = EliminateVarsSimplify[Eqs][{{jsys, jrules}, js}];
            {jsys, jrules} = CleanEqualities[{jsys,jrules}];
            AssociateTo[uvalues, jrules];
        ];
        uvalues
    ]

EliminateVarsStep[Eqs_][{{system_, rules_}, us_, persistus_}] :=
    Module[ {
        subsys,
        newsys,
        newrules
        },
        subsys = Select[system, Function[x, !(FreeQ[x, #])]]& /@ us;
        subsys = And @@ subsys;
        subsys = BooleanConvert[Reduce[ Reduce @ #, Reals]]&[subsys];
        {subsys, newrules} = CleanEqualitiesOperator[Eqs][{subsys, rules}];
        newsys = (system && subsys) /. newrules;
        (*newsys = Reduce[#, Reals]& /@ newsys;*)
        newsys = (BooleanConvert[newsys, "CNF"]);
        {{newsys, newrules}, us, Intersection[us, getVar[newsys]]}
    ]
    
EliminateVarsStep[Eqs_][{{system_, rules_}, {}, {}}] :=
    {{system, rules}, {}, {}}
 
EliminateVarsStep[Eqs_][{{True, rules_}, us_, persistus_}] :=
    {{True, rules}, us, persistus}
 
EliminateVars[Eqs_][{{system_, rules_}, us_}] :=
    EliminateVars[Eqs][{{system /. rules, rules}, us, {}}]     

EliminateVars[Eqs_][{{system_, rules_}, us_, persistus_}] :=
    FixedPoint[EliminateVarsStep[Eqs], {{system /. rules, rules}, us, persistus}]     

FixedReduce3[Eqs_Association][<||>] :=
    CriticalCongestionSolver2[Eqs]   

FixedReduce3[Eqs_Association][{}] :=
    CriticalCongestionSolver2[Eqs]   
    
FixedReduce3[Eqs_Association][approxrules_Association] :=
    Module[ {(*auxsys, auxsol, error, TOL = 10^-10,*) system, 
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
        (*nonlinear = And @@ (MapThread[(Equal[#1,#2]) &, {Eqs["Nlhs"], RoundValues[Eqs["Nrhs"]/. approxrules]}]);
        nonlinear = nonlinear /. Eqs["InitRules"];*)
        (*Lets use "RuleNonCritical" instead.*)
(*        Print[Eqs["Nrhs"]];
        Print["nonlinear ewqs\n",nonlinear];
        Print["replacement rules: \n", Eqs["costpluscurrents"]/.approxrules];*)
        AllIneqs = Lookup[Eqs, "AllIneq", Print["No inequalities to solve."];
                                          Return["a"]];
        AllOrs = Lookup[Eqs, "AllOr", Print["No alternatives to solve. \n",AllOrs];
                                      Return["b"]];
        jrules = js/.approxrules;
        Print[Eqs["EqNonCritical"]/.Eqs["costpluscurrents"]/.jrules];
        Print["non critical rules with j values:\n",rules];
        (*Replace critical congestion rules!*)
        AllOrs = AllOrs /. rules;
        AllIneqs = AllIneqs /. rules;
        
        (*try to simplify bit by bit:*)
        system = AllIneqs && AllOrs;
        rvars = Variables[Join[us,js]/.rules];
        Print["EliminateVarsSimplify for the us"];
        {{system, rules}, aux, newus} = EliminateVarsSimplify[Eqs][{{system, rules}, (*new*)us, {}}];
        rules = Expand/@rules;
        If[ Variables[us/.rules] === {},
            Print["Finished with the us!"],
            Print["Not finished with the us!"];
            {{system, rules}, aux, newus} = EliminateVars[Eqs][{{system, rules}, (*new*)us, {}}];
            Print["Reducing further (NewReduce)..."(*, system*)];
            system = NewReduce[system];
            Print["Simplifying..."];
            system = Simplify @ system;
            {system, rules} = CleanEqualities[{system, rules}];
            uvalues = AssociationThread[us, us /. rules];
            jrules = AssociationThread[js, js /. rules];
            AssociateTo[uvalues, jrules];
            If[ Variables[us/.uvalues] =!= {},
                Print["(Maybe the) system is not feasible!\n
            \tCheck the paper: The current method for stationary mean-field games on networks"]
            ];
        ];
        If[ Variables[js/.rules] === {},
            Print["Done!\n"];
            uvalues = AssociationThread[Join[us, js], Join[us, js] /. rules],
            Print["Finish with the js"];
            uvalues = AssociationThread[us, us /. rules];
            (*TODO Fix this! Use rules or correct the equations! not the criticalcase!*)
            jsys = (*(Eqs["EqCriticalCase"] /. uvalues) &&*) Eqs["EqPosJs"] && Eqs["EqCurrentCompCon"];
            (**Retrieve js values already defined*)
            (*Keeping the order for the js Association*)
            jrules = AssociationThread[js, js/.Select[AssociationThread[js, js/.rules], NumericQ]];
            Print["jrules: \n", jrules];
            {jsys, jrules} = CleanEqualities[{jsys, jrules}];
            {{jsys, jrules}, aux, newjs} = EliminateVarsSimplify[Eqs][{{jsys, jrules}, js}];
            {jsys, jrules} = CleanEqualities[{jsys,jrules}];
            AssociateTo[uvalues, jrules];
        ];
        Print["returned uvalues: \n",Values@uvalues//N];
        uvalues
    ]
 (*    Pause[30];       

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
        (*auxsol*)
 *)
 
 

End[]