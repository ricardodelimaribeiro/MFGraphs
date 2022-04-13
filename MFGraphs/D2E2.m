(* Wolfram Language package *)


ConsistentSwithingCosts::usage = 
"ConsistentSwithingCosts[switching costs][one switching cost]  returns the condition for this switching cost to satisfy the triangle inequality"
ConsistentSwithingCosts[sc_][{a_, b_, c_, S_}] :=
    Module[ {origin, destination, bounds},
        origin = Cases[sc, {a, b, _, _}];
        destination = Cases[sc, {_, b, c, _}];
        bounds = Outer[Plus, Last /@ origin, Last /@ destination] // Flatten;
        And @@ (#>=S &) /@ bounds
    ];
    
IsSwitchingCostConsistent::usage = 
"IsSwitchingCostConsistent[list of switching costs] is True if all switching costs satisfy the triangle inequality"
IsSwitchingCostConsistent[SC_] :=
    Simplify[ConsistentSwithingCosts[SC] /@ SC]
    
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

Data2Equations::usage =
"Data2Equations[Data] returns the equations, inequalities, and alternatives associated to the Data"
Data2Equations[Data_Association] :=
    Module[ { VL = Lookup[Data, "Vertices List", {}]
        , AM = Lookup[Data, "Adjacency Matrix", {}]
        , EVC = Lookup[Data, "Entrance Vertices and Currents", {}]
        , EVTC = Lookup[Data, "Exit Vertices and Terminal Costs", {}]
        , SC = Lookup[Data, "Switching Costs", {}]
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
        , costpluscurrents
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
        
        },
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
        costpluscurrents = Table[Symbol["cpc" <> ToString[k]], {k, 1, Length @ BEL}];
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
        
        
        ModuleVars = (*list of all module variables, except for ModuleVars*)ModuleVars;
        Join[Data, AssociationThread[ToString /@ ModuleVars, ModuleVars]]
    ];
    
GetKirchhoffMatrix[Eqs_] :=
    Module[ { Kirchhoff
        , EqEntryIn = Lookup[Eqs, "EqEntryIn", True]
        , BalanceGatheringCurrents = Lookup[Eqs, "BalanceGatheringCurrents", {}]
        , BalanceSplittingCurrents = Lookup[Eqs, "BalanceSplittingCurrents", {}]
        },
        Kirchhoff = Join[EqEntryIn, (# == 0 & /@ (BalanceGatheringCurrents + BalanceSplittingCurrents))];
    ];    
		