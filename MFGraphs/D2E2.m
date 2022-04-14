(* Wolfram Language package *)


ConsistentSwithingCosts::usage = 
"ConsistentSwithingCosts[switching costs][one switching cost]  returns the condition for this switching cost to satisfy the triangle inequality"
ConsistentSwithingCosts[sc_][{a_, b_, c_, S_}] :=
    Module[ {origin, destination, bounds},
        origin = Cases[sc, {a, b, _, _}];
        destination = Cases[sc, {_, b, c, _}];
        bounds = Outer[Plus, Last /@ origin, Last /@ destination] // Flatten;
        Print[bounds];
        And @@ (#>=S &) /@ bounds
    ];
    
IsSwitchingCostConsistent::usage = 
"IsSwitchingCostConsistent[list of switching costs] is True if all switching costs satisfy the triangle inequality"
IsSwitchingCostConsistent[SC_] :=
    And @@ Simplify[ConsistentSwithingCosts[SC] /@ SC]
    
AtHead[DirectedEdge[a_, b_]] :=
    {b, DirectedEdge[a, b]}

AtTail[DirectedEdge[a_, b_]] :=
    {a, DirectedEdge[a, b]}

TransitionsAt[G_, k_] :=
    Prepend[#, k] & /@ Permutations[IncidenceList[G, k], {2}]

RoundValues[x_?NumberQ] :=
    Round[x, 10^-10]
    
RoundValues[Rule[a_, b_]] :=
    Rule[a, RoundValues[b]]
    
RoundValues[x_List] :=
    RoundValues /@ x
    
RoundValues[x_Association] :=
    RoundValues /@ x

triple2path[{a_, b_, c_}, G_] :=
    Module[ {EL = EdgeList[G]},
        If[ SubsetQ[AdjacencyList[G, b], {a, c}],
            If[ MemberQ[EL, DirectedEdge[a, b]],
                If[ MemberQ[EL, DirectedEdge[b, c]],
                    {b, DirectedEdge[a, b] , DirectedEdge[b, c]},
                    {b, DirectedEdge[a, b], DirectedEdge[c, b]}
                ],
                If[ MemberQ[EL, DirectedEdge[b, a]],
                    If[ MemberQ[EL, DirectedEdge[b, c]],
                        {b, DirectedEdge[b, a], DirectedEdge[b, c]},
                        {b, DirectedEdge[b, a], DirectedEdge[c, b]}
                    ]
                ]
            ],
            StringJoin["\nThere is no path from ", ToString[a], " to ", ToString[b], " to ", ToString[c]]
        ]
    ]
  
path2triple[{a_, b_ \[DirectedEdge] c_, d_ \[DirectedEdge] e_} -> S_] :=
    {Select[{b, c}, # =!= a &], a, Select[{d, e}, # =!= a &], S} // Flatten;
       
CurrentCompCon[jvars_][a_ \[DirectedEdge] b_] :=
    jvars[{a, a \[DirectedEdge] b}] == 0 || jvars[{b, a \[DirectedEdge] b}] == 0;

CurrentSplitting[AllTransitions_][{c_, a_ \[DirectedEdge] b_}] :=
    Select[AllTransitions, (Take[#, 2] == {c, a \[DirectedEdge] b}) &];

CurrentGathering[AllTransitions_][{c_, a_ \[DirectedEdge] b_}] :=
    Select[AllTransitions, (Part[#, {1, 3}] == OtherWay[{c, a \[DirectedEdge] b}]) &];

TransitionCompCon[jtvars_][{v_, edge1_, edge2_}] :=
    jtvars[{v, edge1, edge2}] == 0 || jtvars[{v, edge2, edge1}] == 0;
        
IncomingEdges[FG_][k_] :=
    {k, #1} & /@ IncidenceList[FG,k];

OtherWay[{c_, DirectedEdge[a_, b_]}] :=
    {If[ c === a,
         b,
         a
     ], DirectedEdge[a, b]}

OutgoingEdges[FG_][k_] :=
    OtherWay /@ ({k, #} & /@ IncidenceList[FG, k]);

ExitRules[uvars_,ExitCosts_][a_ \[DirectedEdge] b_] :=
    Total[uvars /@ {{b, DirectedEdge[a, b]}}] -> ExitCosts[b];

ExitCurrents[jvars_][a_ \[DirectedEdge] b_] :=
    jvars @ {a, DirectedEdge[a, b]} -> 0;

Transu[uvars_,SwitchingCosts_][{v_, edge1_, edge2_}] :=
    uvars[{v,edge1}] <= uvars[{v,edge2}] + SwitchingCosts[{v,edge1,edge2}];
     
Compu[jtvars_,uvars_,SwitchingCosts_][{v_, edge1_, edge2_}] :=
    (jtvars[{v, edge1, edge2}] == 0) || uvars[{v,edge2}]-uvars[{v,edge1}] + SwitchingCosts[{v,edge1,edge2}] == 0;
        
Data2Equations::usage =
"Data2Equations[Data] returns the equations, inequalities, and alternatives associated to the Data"
Data2Equations[Data_Association] :=
    Module[ { VL
        , AM
        , EVC
        , EVTC
        , SC
        , BG
        , EntranceVertices
        , InwardVertices
        , ExitVertices
        , OutwardVertices
        , InEdges
        , OutEdges
        , AuxiliaryGraph
        , FG
        , EL
        , BEL
        , FVL
        , jargs    
        , uargs
        , AllTransitions
        , EntryArgs
        , EntryDataAssociation
        , ExitCosts
        , js
        , jvars
        , us
        , uvars
        , jts
        , jtvars
        , SignedCurrents
        , SwitchingCosts
        , EqPosJs
        , EqPosJts
        , EqCurrentCompCon
        , EqTransitionCompCon
        , NoDeadEnds
        , EqBalanceSplittingCurrents
        , BalanceSplittingCurrents
        , NoDeadStarts
        , RuleBalanceGatheringCurrents
        , BalanceGatheringCurrents
        , EqBalanceGatheringCurrents
        , EqEntryIn
        , RuleEntryOut
            , RuleExitCurrentsIn
            , RuleExitValues
            , EqValueAuxiliaryEdges
            , OutRules
            , InRules
            , EqSwitchingByVertex
            , EqCompCon
            , Nlhs
            , ModuleVars
            , ModuleVarsNames
        },
        VL = Lookup[Data, "Vertices List", {}];
        AM = Lookup[Data, "Adjacency Matrix", {}];
        EVC = Lookup[Data, "Entrance Vertices and Currents", {}];
        EVTC = Lookup[Data, "Exit Vertices and Terminal Costs", {}];
        SC = Lookup[Data, "Switching Costs", {}];
            (****Graph stuff****)
        BG = AdjacencyGraph[VL,AM, VertexLabels -> "Name", DirectedEdges -> True];
        EntranceVertices = First /@ EVC;
        ExitVertices = First /@ EVTC;
        Clear["en*","ex*"];
        
        (*InwardVertices defines auxiliary vertices for the entrance vertices*)
        InwardVertices = AssociationThread[EntranceVertices, Symbol["en" <> ToString[#]] & /@ EntranceVertices];
        OutwardVertices = AssociationThread[ExitVertices, Symbol["ex" <> ToString[#]] & /@ ExitVertices];
        
        (*InEdges defines auxiliary arguments for the entrance vertices*)
        InEdges = MapThread[DirectedEdge, {InwardVertices /@ EntranceVertices, EntranceVertices}];
        OutEdges = MapThread[DirectedEdge, {ExitVertices, OutwardVertices /@ ExitVertices}];
        AuxiliaryGraph = Graph[Join[InEdges, OutEdges], VertexLabels -> "Name", GraphLayout -> "SpringEmbedding"];
        FG = EdgeAdd[BG, Join[InEdges, OutEdges]];
        EL = EdgeList[FG];
        BEL = EdgeList[BG];
        FVL = VertexList[FG];
        
        (*arguments*)
        jargs = Flatten[#, 1] & @ ({AtTail@#, AtHead@#} & /@ EL);
        uargs = jargs;
        AllTransitions = TransitionsAt[FG, #] & /@ FVL // Catenate(*at vertex from first edge to second edge*);
        EntryArgs = AtHead /@ ((EdgeList[AuxiliaryGraph, DirectedEdge[_,#] ] & /@ (First /@ EVC)) // Flatten[#, 1] &);
        EntryDataAssociation = RoundValues @ AssociationThread[EntryArgs, Last /@ EVC];
        ExitCosts = AssociationThread[OutwardVertices /@ (First /@ EVTC), Last /@ EVTC];
        
        (*variables*)
        js = Table[Symbol["j" <> ToString[k]], {k, 1, Length @ jargs}];
        jvars = AssociationThread[jargs, js];
        jts = Table[Symbol["jt" <> ToString[k]], {k, 1, Length @ AllTransitions}];
        jtvars = AssociationThread[AllTransitions, jts];
        us = Table[Symbol["u" <> ToString[k]], {k, 1, Length @ uargs}];
        uvars = AssociationThread[uargs, us];
        SignedCurrents =   AssociationThread[BEL, (jvars[AtHead[#]] - jvars[AtTail[#]] &) /@ BEL];
        Print["D2E: Variables are all set"];
        
        (*Elements of the system*)
        (*Swithing cost is initialized with 0. AssociationThread associates the last association!*)
        SwitchingCosts = AssociationThread[Join[AllTransitions, triple2path[Take[#, 3], FG] & /@ SC], Join[0&/@ AllTransitions, Last[#] & /@ SC]];
        EqPosJs = And @@ ( #>=0& /@ Join[jvars]);(*Inequality*)
        EqPosJts = And @@ ( #>=0& /@ Join[jtvars]);(*Inequality*)
        EqCurrentCompCon = And @@ (CurrentCompCon[jvars] /@ EL);(*Or*)
        EqTransitionCompCon = And @@ ((Sort /@ TransitionCompCon[jtvars] /@ AllTransitions) // Union);(*Or*)
        
        (*Balance Splitting Currents in the full graph*)
        NoDeadEnds = IncomingEdges[FG] /@ VL // Flatten[#, 1] &;
        BalanceSplittingCurrents = ((jvars[#] - Total[jtvars /@ CurrentSplitting[AllTransitions][#]]) & /@ NoDeadEnds);
        EqBalanceSplittingCurrents = Simplify /@ (And @@ ((#==0)& /@ BalanceSplittingCurrents));(*Equal*)
        
        (*Gathering currents in the inside of the basic graph*)
        NoDeadStarts = OutgoingEdges[FG] /@ VL // Flatten[#, 1] &;
        RuleBalanceGatheringCurrents = (jvars[#] -> Total[jtvars /@ CurrentGathering[AllTransitions][#]]) & /@ NoDeadStarts;(*Rule*)
        
        (*get equations for the exit currents at the entry vertices*)
        BalanceGatheringCurrents = ((-jvars[#] + Total[jtvars /@ CurrentGathering[AllTransitions][#]]) & /@ NoDeadStarts);
        EqBalanceGatheringCurrents = Simplify/@(And @@ (#==0&/@ BalanceGatheringCurrents));
        
        (*Incoming currents*)
        EqEntryIn = (jvars[#] == EntryDataAssociation[#]) & /@ (AtHead /@ InEdges);(*List of Equals*)
        
        (*Outgoing currents at entrances*)
        RuleEntryOut = (jvars[#] -> 0) & /@ (AtTail /@ InEdges);(*Rule*)
        RuleExitCurrentsIn = ExitCurrents[jvars]/@ OutEdges;(*Rule*)
        
        (*Exit values at exit vertices*)
        RuleExitValues = ExitRules[uvars,ExitCosts] /@ OutEdges;(*Rule*)
        
        (*The value function on the auxiliary edges is constant and equal to the exit cost.*)
        EqValueAuxiliaryEdges = And @@ ((uvars[AtTail[#]] == uvars[AtHead[#]]) & /@ Join[InEdges, OutEdges]);(*Equal*)(*use ToRules to get the rules*)
        OutRules = Rule[#, Infinity] & /@ (Outer[Flatten[{AtTail[#1], #2}] &, OutEdges, EL] // Flatten[#, 1] &);
        InRules = Rule[#, Infinity] & /@ (Outer[{#2[[2]], #1, #2} &, IncidenceList[FG, #] & /@ EntranceVertices, InEdges] // Flatten[#, 2] &);
        
        (*Switching condition equations*)
        EqSwitchingByVertex = Transu[uvars,SwitchingCosts]/@TransitionsAt[FG, #] & /@ VL;
        EqCompCon = And @@ Compu[jtvars,uvars,SwitchingCosts] /@ AllTransitions;(*Or*)
        Nlhs = Flatten[uvars[AtHead[#]] - uvars[AtTail[#]] + SignedCurrents[#] & /@ BEL];
        ModuleVars = (*list of all module variables, except for ModuleVars*)
        { VL
        , AM
        , EVC
        , EVTC
        , SC
        , BG
        , EntranceVertices
        , InwardVertices
        , ExitVertices
        , OutwardVertices
        , InEdges
        , OutEdges
        , AuxiliaryGraph
        , FG
        , EL
        , BEL
        , FVL
        , jargs    
        , uargs
        , AllTransitions
        , EntryArgs
        , EntryDataAssociation
        , ExitCosts
        , js
        , jvars
        , us
        , uvars
        , jts
        , jtvars
        , SignedCurrents
        , SwitchingCosts
        , EqPosJs
        , EqPosJts
        , EqCurrentCompCon
        , EqTransitionCompCon
        , NoDeadEnds
        , EqBalanceSplittingCurrents
        , BalanceSplittingCurrents
        , NoDeadStarts
        , RuleBalanceGatheringCurrents
        , BalanceGatheringCurrents
        , EqBalanceGatheringCurrents
        , EqEntryIn
        , RuleEntryOut
            , RuleExitCurrentsIn
            , RuleExitValues
            , EqValueAuxiliaryEdges
            , OutRules
            , InRules
            , EqSwitchingByVertex
            , EqCompCon
            , Nlhs
            };
        ModuleVarsNames = {"VL", "AM", "EVC", "EVTC", "SC", "BG", "EntranceVertices", \
"InwardVertices", "ExitVertices", "OutwardVertices", "InEdges", \
"OutEdges", "AuxiliaryGraph", "FG", "EL", "BEL", "FVL", "jargs", \
"uargs", "AllTransitions", "EntryArgs", "EntryDataAssociation", \
"ExitCosts", "js", "jvars", "us", "uvars", "jts", "jtvars", \
        "SignedCurrents", "SwitchingCosts", "EqPosJs", \
"EqPosJts", "EqCurrentCompCon", "EqTransitionCompCon", "NoDeadEnds", \
"EqBalanceSplittingCurrents", "BalanceSplittingCurrents", \
"NoDeadStarts", "RuleBalanceGatheringCurrents", \
"BalanceGatheringCurrents", "EqBalanceGatheringCurrents", \
"EqEntryIn", "RuleEntryOut", "RuleExitCurrentsIn", "RuleExitValues", \
"EqValueAuxiliaryEdges", "OutRules", "InRules", \
"EqSwitchingByVertex", "EqCompCon", "Nlhs"};
        Print[ModuleVarsNames];
        Join[Data, AssociationThread[ModuleVarsNames, ModuleVars]]
    ];
    
GetKirchhoffMatrix[Eqs_] :=
    Module[ { Kirchhoff
        , EqEntryIn = Lookup[Eqs, "EqEntryIn", True]
        , BalanceGatheringCurrents = Lookup[Eqs, "BalanceGatheringCurrents", {}]
        , BalanceSplittingCurrents = Lookup[Eqs, "BalanceSplittingCurrents", {}]
        , RuleExitCurrentsIn = Lookup[Eqs, "RuleExitCurrentsIn", {}]
        , RuleEntryOut = Lookup[Eqs, "RuleEntryOut", {}]
        , BM
        , KM
        , vars
        },
        Kirchhoff = Join[EqEntryIn, (# == 0 & /@ (BalanceGatheringCurrents + BalanceSplittingCurrents))];
        Kirchhoff = Kirchhoff /. Join[RuleExitCurrentsIn,RuleEntryOut];
        {BM,KM} = CoefficientArrays[Kirchhoff, vars = Variables[Kirchhoff /. Equal -> Plus]];
        Print["The matrices B and K are: \n",MatrixForm/@{-BM,KM},"\nThe order of the variables is \n", vars];
        {-BM,KM}
    ];    
    
MFGPreprocessing[Eqs_] :=
    Module[ {InitRules
        , RuleBalanceGatheringCurrents = Lookup[Eqs, "RuleBalanceGatheringCurrents", {}]
        , RuleEntryIn = Lookup[Eqs, "RuleEntryIn", True]
        , RuleEntryOut = Lookup[Eqs, "RuleEntryOut", True]
        , TrueEq
        , EqBalanceGatheringCurrents = Lookup[Eqs, "EqBalanceGatheringCurrents", True]
        , RuleExitCurrentsIn = Lookup[Eqs, "RuleExitCurrentsIn", True]
        , RuleExitValues = Lookup[Eqs, "RuleExitValues", True]
        , EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", True]
        , EqSwitchingConditions
        , EqSwitchingByVertex = Lookup[Eqs, "EqSwitchingByVertex", True]
        , EqCompCon = Lookup[Eqs, "EqCompCon", True]
        , EqBalanceSplittingCurrents = Lookup[Eqs, "EqBalanceSplittingCurrents", True]
        , AllOr
        , EqCurrentCompCon = Lookup[Eqs, "EqCurrentCompCon", True]
        , EqTransitionCompCon = Lookup[Eqs, "EqTransitionCompCon", True]
        , EqPosCon = Lookup[Eqs, "EqPosCon", True]
        , AllIneq
        , newStuff
    },
    (*Check switching costs for consistency*)
    
    (*First rules: these have some j in terms of jts*)
        InitRules = Association[RuleBalanceGatheringCurrents];
        AssociateTo[InitRules, Join[RuleEntryIn, RuleEntryOut]];
        (*I think this does nothing!!!
        Include Gathering currents information in the rules*)
        {TrueEq, InitRules} = CleanEqualities[{EqBalanceGatheringCurrents, InitRules}];
        AssociateTo[InitRules, RuleExitCurrentsIn];
        AssociateTo[InitRules, RuleExitValues];
        {TrueEq, InitRules} = CleanEqualities[{EqValueAuxiliaryEdges, InitRules}];
    (**)
    		
    		{EqSwitchingConditions, InitRules} = CleanEqualities[{EqSwitchingByVertex, InitRules}];
    		{EqCompCon, InitRules} = CleanEqualities[{EqCompCon, InitRules}];
    		{TrueEq, InitRules} = CleanEqualities[{EqBalanceSplittingCurrents, InitRules}];
    		
    		AllOr = EqCurrentCompCon && EqTransitionCompCon && EqCompCon;
        AllOr = AllOr/.InitRules;
        AllOr = BooleanConvert[Simplify /@ AllOr, "CNF"];
        
        {AllOr, InitRules} = CleanEqualities[{AllOr, InitRules}];
        
        EqPosCon = EqPosCon /. InitRules;
        
        EqSwitchingConditions = EqSwitchingConditions/.InitRules;
        EqSwitchingConditions = Simplify[EqSwitchingConditions];
        AllOr = AllOr && Select[EqSwitchingConditions,Head[#] === Or &];
        AllOr = BooleanConvert[AllOr,"CNF"];
        AllIneq = EqPosCon && Select[EqSwitchingConditions,Head[#] =!= Or &];
        AllIneq = Simplify[AllIneq];
    
    (*return something!!!*)
    Join[Eqs, newStuff]
    
        ];
CriticalCongestionSolver[Eqs_] :=
    Module[ {PreEqs
    },
        PreEqs = MFGPreprocessing[Eqs];
        (*Solve the updated (with rules from preprocessing) edge equations*)
        (*Replace in preprocessed system and rules and possibly use cleanequalities or newreduce...*)
        Solve[PreEqs[""]]
    ];
    
(*EqEliminator[system_, rules_List] :=
    EqEliminator[system, Association[rules]]
*)
EqEliminator[{sys_, rules_List}] :=
	EqEliminator[{sys, Association @ rules}]
	
EqEliminator[{sys_, rules_Association}] :=
    Module[ {EE = True, ON = True, newrules,groups, system = sys /. rules},
        Which[
        Head[system] === And,    
        		(*separate equalities from the rest*)
        		groups = GroupBy[List@@system, Head[#]===Equal&];
        		EE = And@@Lookup[groups, True, Return[{system, rules},Module]];
        		ON = And@@Lookup[groups, False, True],
        Head[system] === Equal,
        		EE = system,
        True,(*Or or Inequality*)
        		ON = Simplify @ system
        ];
        newrules = First @ Solve[EE] // Quiet; (*The reason we use Quiet is:  Solve::svars: Equations may not give solutions for all "solve" variables.*)
        newrules = Join[rules /. newrules, Association @ newrules];
        {ON , Expand /@ newrules}
    ]

CleanEqualities[system_, rules_] :=
    FixedPoint[EqEliminator, {system, rules}]

CleanEqualities[system_List, rules_] :=
    Module[ {systems, ruless, aux},
        aux = CleanEqualities[#, rules] &/@ system;
        systems = First /@ aux;
        systems = And @@ systems;
        ruless = Last /@ aux;
        ruless = Join @@ ruless;
        systems = systems /. ruless;
        ruless = ruless /. ruless;
        {systems, ruless}
    ]
	
