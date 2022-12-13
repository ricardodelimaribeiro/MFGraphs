(*Wolfram Language package*)
ConsistentSwithingCosts::usage = "ConsistentSwithingCosts[switchingcosts][{a,b,c}->S]  
returns True if S, the cost of switching from the edge ab to cb, is smaller than any other combination, such as, ab to bd and then from bd to dc.
returns the condition for this switching cost to satisfy the triangle inequality when S, and the other switching costs too, does not have a numerical value.";

ConsistentSwithingCosts[sc_][{a_, b_, c_} -> S_] :=
    Module[ {origin, destination, bounds},
        origin = Cases[sc, HoldPattern[{a, b, _} -> _]];
        origin = DeleteCases[origin, {a, b, c} -> S];
        If[ origin =!= {},
            destination = {a, Part[First[#], 3], c} & /@ origin;
            bounds = ((S <= Last[#] + Association[sc][{a, Part[First[#], 3], c}]) & /@ origin);
            And @@ bounds,
            True
        ]
    ];

IsSwitchingCostConsistent::usage = "IsSwitchingCostConsistent[List of \
switching costs] is True if all switching costs satisfy the triangle \
inequality. If some switching costs are symbolic, then it returns the consistency conditions."
IsSwitchingCostConsistent[SC_] :=
    And @@ Simplify[ConsistentSwithingCosts[SC] /@ SC]

AtHead::usage = 
"AtHead[DirectedEdge[a,b] returns {b,DirectedEdge[a,b]}. This is a way o selecting the edge with the orientation. "

AtHead[DirectedEdge[a_, b_]] :=
    {b, DirectedEdge[a, b]}

AtTail::usage = 
"AtTail[DirectedEdge[a,b] returns {a,DirectedEdge[a,b]}. This is a way o selecting the edge with the orientation. "

AtTail[DirectedEdge[a_, b_]] :=
    {a, DirectedEdge[a, b]}

TransitionsAt[G_, k_] :=
    Prepend[#, k] & /@ Permutations[IncidenceList[G, k], {2}]

triple2path[{a_, b_, c_}, G_] :=
    Module[ {EL = EdgeList[G]},
        If[ SubsetQ[AdjacencyList[G, b], {a, c}],
            If[ MemberQ[EL, DirectedEdge[a, b]],
                If[ MemberQ[EL, DirectedEdge[b, c]],
                    {b, DirectedEdge[a, b], 
                    DirectedEdge[b, c]},
                    {b, DirectedEdge[a, b], 
                    DirectedEdge[c, b]}
                ],
                If[ MemberQ[EL, DirectedEdge[b, a]],
                    If[ MemberQ[EL, DirectedEdge[b, c]],
                        {b, DirectedEdge[b, a], 
                        DirectedEdge[b, c]},
                        {b, DirectedEdge[b, a], 
                        DirectedEdge[c, b]}
                    ]
                ]
            ],
            StringJoin["\nThere is no path from ", ToString[a], " to ", 
             ToString[b], " to ", ToString[c]]
        ]
    ]

path2triple[{a_, b_ \[DirectedEdge] c_, d_ \[DirectedEdge] e_} -> 
    S_] :=
    {Select[{b, c}, # =!= a &], a, Select[{d, e}, # =!= a &], 
    S} // Flatten;

CurrentCompCon[jvars_][a_ \[DirectedEdge] b_] :=
    jvars[{a, a \[DirectedEdge] b}] == 0 || 
     jvars[{b, a \[DirectedEdge] b}] == 0;

CurrentSplitting[AllTransitions_][{c_, a_ \[DirectedEdge] b_}] :=
    Select[AllTransitions, (Take[#, 2] == {c, a \[DirectedEdge] b}) &];

CurrentGathering[AllTransitions_][{c_, a_ \[DirectedEdge] b_}] :=
    Select[
     AllTransitions, (Part[#, {1, 3}] == 
        OtherWay[{c, a \[DirectedEdge] b}]) &];

TransitionCompCon[jtvars_][{v_, edge1_, edge2_}] :=
    jtvars[{v, edge1, edge2}] == 0 || jtvars[{v, edge2, edge1}] == 0;

IncomingEdges[FG_][k_] :=
    {k, #1} & /@ IncidenceList[FG, k];

OtherWay[{c_, DirectedEdge[a_, b_]}] :=
    {If[ c === a,
         b,
         a
     ], 
    DirectedEdge[a, b]}

OutgoingEdges[FG_][k_] :=
    OtherWay /@ ({k, #} & /@ IncidenceList[FG, k]);

BackTransition[{k_, edge1_, edge2_}] :=
    {k, edge2, edge1}

ExitRules[uvars_, ExitCosts_][a_ \[DirectedEdge] b_] :=
    Total[uvars /@ {{b, DirectedEdge[a, b]}}] -> ExitCosts[b];

ExitCurrents[jvars_][a_ \[DirectedEdge] b_] :=
    jvars@{a, DirectedEdge[a, b]} -> 0;
Transu::usage = 
"Transu[u, SC][{v,e1,e2}] returns the optimality condition at the vertex v related to switching form e1 to e2. Namely, 
u[{v, e1}] <= u[{v, e2}] + SC[{v, e1, e2}]"
Transu[uvars_, SwitchingCosts_][{v_, edge1_, edge2_}] :=
    uvars[{v, edge1}] <= 
     uvars[{v, edge2}] + SwitchingCosts[{v, edge1, edge2}];

Compu::usage =
"Compu[jt,u,SC][{v,e1,e2}] returns the complementarity condition:
(jt[{v, e1, e2}] == 0) || (u[{v, e2}] - u[{v, e1}] + SC[{v, e1, e2}] == 0)"
Compu[jtvars_, uvars_, SwitchingCosts_][{v_, edge1_, 
    edge2_}] :=
    (jtvars[{v, edge1, edge2}] == 0) || 
    (uvars[{v, edge2}] - uvars[{v, edge1}] + SwitchingCosts[{v, edge1, edge2}] == 0);

Data2Equations::usage = "Data2Equations[Data] returns the equations, \
inequalities, and alternatives associated to the Data. "
Data2Equations[Data_Association] :=
    Module[ {VL, AM, EVC, EVTC, SC, BG, EntranceVertices, InwardVertices,
      ExitVertices, OutwardVertices, InEdges, OutEdges, AuxiliaryGraph,
      FG, EL, BEL, FVL, jargs, uargs, AllTransitions, EntryArgs, 
      EntryDataAssociation, ExitCosts, js, jvars, us, uvars, jts, 
      jtvars, SignedCurrents, SwitchingCosts, EqPosJs, EqPosJts, 
      EqCurrentCompCon, EqTransitionCompCon, NoDeadEnds, 
      EqBalanceSplittingCurrents, BalanceSplittingCurrents, 
      NoDeadStarts, RuleBalanceGatheringCurrents, 
      BalanceGatheringCurrents, EqBalanceGatheringCurrents, EqEntryIn, 
      RuleEntryOut, RuleExitCurrentsIn, RuleExitValues, 
      EqValueAuxiliaryEdges, OutRules, InRules, EqSwitchingByVertex, 
      EqCompCon, Nlhs, ModuleVars, ModuleVarsNames, LargeCases, 
      LargeSwitchingTransitions, ZeroRun, CostArgs, Nrhs, 
      consistentCosts, costpluscurrents, EqGeneral, CostArgs2},
        VL = Lookup[Data, "Vertices List", {}];
        AM = Lookup[Data, "Adjacency Matrix", {}];
        EVC = Lookup[Data, "Entrance Vertices and Currents", {}];
        EVTC = Lookup[Data, "Exit Vertices and Terminal Costs", {}];
        SC = Lookup[Data, "Switching Costs", {}];
        
        (****Graph stuff****)
        BG = AdjacencyGraph[VL, AM, VertexLabels -> "Name", DirectedEdges -> True];
        EntranceVertices = First /@ EVC;
        ExitVertices = First /@ EVTC;
        Clear["en*", "ex*"];
        (*InwardVertices defines auxiliary vertices for the entrance \
     		vertices*)
        InwardVertices = AssociationMap[Symbol["en" <> ToString[#]] &, EntranceVertices];
        OutwardVertices = AssociationMap[Symbol["ex" <> ToString[#]] &, ExitVertices];
        
        (*InEdges defines auxiliary arguments for the entrance vertices*)
        InEdges = MapThread[DirectedEdge, {InwardVertices /@ EntranceVertices, EntranceVertices}];
        OutEdges = MapThread[DirectedEdge, {ExitVertices, OutwardVertices /@ ExitVertices}];
        AuxiliaryGraph = Graph[Join[InEdges, OutEdges], VertexLabels -> "Name", GraphLayout -> "SpringEmbedding"];
        FG = EdgeAdd[BG, Join[InEdges, OutEdges]];
        EL = EdgeList[FG];
        BEL = EdgeList[BG];
        FVL = VertexList[FG];
        
        (*arguments*)
        jargs = Flatten[#, 1] &@({AtTail@#, AtHead@#} & /@ EL);
        uargs = jargs;
        AllTransitions = TransitionsAt[FG, #] & /@ FVL // Catenate(*at vertex from first edge to second edge*);
        EntryArgs = AtHead /@ ((EdgeList[AuxiliaryGraph, DirectedEdge[_, #]] & /@ (First /@ EVC)) // Flatten[#, 1] &);
        ZeroRun = AssociationMap[0 &][AtHead /@ Join[InEdges, OutEdges]];
        EntryDataAssociation = RoundValues@AssociationThread[EntryArgs, Last /@ EVC];
        ExitCosts = AssociationThread[OutwardVertices /@ (First /@ EVTC), Last /@ EVTC];
        
        (*variables*)
        js = Table[Symbol["j" <> ToString[k]], {k, 1, Length@jargs}];
        jvars = AssociationThread[jargs, js];
        jts = Table[Symbol["jt" <> ToString[k]], {k, 1, Length@AllTransitions}];
        jtvars = AssociationThread[AllTransitions, jts];
        us = Table[Symbol["u" <> ToString[k]], {k, 1, Length@uargs}];
        uvars = AssociationThread[uargs, us];
        SignedCurrents = AssociationMap[jvars[AtHead[#]] - jvars[AtTail[#]] &, EL];
        SwitchingCosts = AssociationMap[0 &, AllTransitions];
        AssociateTo[SwitchingCosts, AssociationThread[triple2path[Take[#, 3], FG] & /@ SC, Last[#] & /@ SC]];
        LargeCases = Join[{Last[#], __, #} & /@ InEdges, {First[#], #, __} & /@ OutEdges];
        LargeSwitchingTransitions = Cases[AllTransitions, #] & /@ LargeCases // Flatten[#, 1] &;
        AssociateTo[SwitchingCosts, AssociationMap[Infinity &, LargeSwitchingTransitions]];
        consistentCosts = IsSwitchingCostConsistent[Normal@SwitchingCosts];
        Which[consistentCosts === False, 
            Print["Switching costs are inconsistent!"];
            Return[ $Failed, Module], 
         consistentCosts =!= True, 
         Print["Switching costs conditions are ", consistentCosts]
        ];
        Clear[CostArgs];
        (*Cost 1: Absolute values of currents plus the product of opposing currents*)
        CostArgs = 
         Join[AssociationMap[Function[x,Abs[SignedCurrents[Last[x]]] + jvars[x]jvars[OtherWay[x]]], jargs], 
             ZeroRun, AssociationMap[jtvars[#] jtvars[BackTransition[#]] &, Keys@SwitchingCosts] + (SwitchingCosts /. {Infinity -> 10^6})];
        (*Cost 1.5 with Hamiltonian: Absolute values of currents plus the product of opposing currents*)
        CostArgs = 
         Join[AssociationMap[Function[x,Abs[Cost[SignedCurrents[Last[x]], Last[x]]] + jvars[x]jvars[OtherWay[x]]], jargs], 
             ZeroRun, AssociationMap[jtvars[#] jtvars[BackTransition[#]] &, Keys@SwitchingCosts] + (SwitchingCosts /. {Infinity -> 10^6})];
        (*Cost 2: sum of absolute values of the opposing currents.*)
        CostArgs2 = 
          Join[AssociationMap[Function[x,10^(-6) + Abs[jvars[x]] + Abs[jvars[OtherWay[x]]]], jargs], 
              ZeroRun, AssociationMap[Abs[jtvars[#]]+ Abs[jtvars[BackTransition[#]]] &, 
        Keys@SwitchingCosts] + (SwitchingCosts /. {Infinity -> 10^6})];
        EqPosJs = And @@ (# >= 0 & /@ Join[jvars]);(*Inequality*)
        EqPosJts = And @@ (# >= 0 & /@ Join[jtvars]);(*Inequality*)
        EqCurrentCompCon = And @@ (CurrentCompCon[jvars] /@ EL);(*Or*)
        EqTransitionCompCon = And @@ ((Sort /@ TransitionCompCon[jtvars] /@ AllTransitions) // Union);(*Or*)(*Balance Splitting Currents in the full graph*)
        NoDeadEnds = IncomingEdges[FG] /@ VL // Flatten[#, 1] &;
        BalanceSplittingCurrents = ((jvars[#] - Total[jtvars /@ CurrentSplitting[AllTransitions][#]]) & /@ NoDeadEnds);
        EqBalanceSplittingCurrents = Simplify /@ (And @@ ((# == 0) & /@ BalanceSplittingCurrents));(*Equal*)(*Gathering currents in \
     		the inside of the basic graph*)
        NoDeadStarts = OutgoingEdges[FG] /@ VL // Flatten[#, 1] &;
        RuleBalanceGatheringCurrents = Association[(jvars[#] -> Total[jtvars /@ CurrentGathering[AllTransitions][#]]) & /@ NoDeadStarts];(*Rule*)(*get equations for the exit currents at \
     		the entry vertices*)
        BalanceGatheringCurrents = ((-jvars[#] + Total[jtvars /@ CurrentGathering[AllTransitions][#]]) & /@ NoDeadStarts);
        EqBalanceGatheringCurrents = Simplify /@ (And @@ (# == 0 & /@ BalanceGatheringCurrents));
        
        (*Incoming currents*)
        EqEntryIn = (jvars[#] == EntryDataAssociation[#]) & /@ (AtHead /@ InEdges);(*List of Equals*)(*Outgoing currents at entrances*)
        RuleEntryOut = Association[(jvars[#] -> 0) & /@ (AtTail /@ InEdges)];(*Rule*)
        RuleExitCurrentsIn = Association[ExitCurrents[jvars] /@ OutEdges];(*Rule*)(*Exit values at exit vertices*)
        RuleExitValues = Association[ExitRules[uvars, ExitCosts] /@ OutEdges];(*Rule*)(*The value function on the auxiliary edges \
     		is constant and equal to the exit cost.*)
        EqValueAuxiliaryEdges = And @@ ((uvars[AtTail[#]] == uvars[AtHead[#]]) & /@ Join[InEdges, OutEdges]);(*Equal*)(*use ToRules to get the rules*)
        OutRules = Rule[#, Infinity] & /@ (Outer[Flatten[{AtTail[#1], #2}] &, OutEdges, 
             EL] // Flatten[#, 1] &);
        InRules = Rule[#, Infinity] & /@ (Outer[{#2[[2]], #1, #2} &, 
             IncidenceList[FG, #] & /@ EntranceVertices, InEdges] // Flatten[#, 2] &);
        EqSwitchingByVertex = Transu[uvars, SwitchingCosts] /@ TransitionsAt[FG, #] & /@ VL;
        EqSwitchingByVertex = (And @@ #) & /@ EqSwitchingByVertex;
        EqSwitchingByVertex = Select[#,FreeQ[Infinity]]&/@EqSwitchingByVertex;
        EqCompCon = And @@ Compu[jtvars, uvars, SwitchingCosts] /@ AllTransitions;(*Or*)
        Nlhs = Flatten[uvars[AtHead[#]] - uvars[AtTail[#]] + SignedCurrents[#] & /@ BEL];
        
        (*SignedCurrents[#] = jvars[AtHead[#]] - jvars[AtTail[#]*)
        Nrhs = Flatten[SignedCurrents[#] - Sign[SignedCurrents[#]] Cost[SignedCurrents[#], #] & /@ BEL];
       	(*stuff to solve the general case faster*)
        costpluscurrents = Table[Symbol["cpc" <> ToString[k]], {k, 1, Length@BEL}];
        EqGeneral = And @@ (MapThread[Equal, {Nlhs, costpluscurrents}]);
        costpluscurrents = AssociationThread[costpluscurrents, Nrhs];
        (*stuff to solve the general case faster*)

        (*list of all module variables, except for ModuleVars*)
        ModuleVars = {VL, AM, EVC, EVTC, SC, BG, 
          EntranceVertices, InwardVertices, ExitVertices, OutwardVertices, 
          InEdges, OutEdges, AuxiliaryGraph, FG, EL, BEL, FVL, jargs, 
          uargs, AllTransitions, EntryArgs, EntryDataAssociation, 
          ExitCosts, js, jvars, us, uvars, jts, jtvars, SignedCurrents, 
          SwitchingCosts, EqPosJs, EqPosJts, EqCurrentCompCon, 
          EqTransitionCompCon, NoDeadEnds, EqBalanceSplittingCurrents, 
          BalanceSplittingCurrents, NoDeadStarts, 
          RuleBalanceGatheringCurrents, BalanceGatheringCurrents, 
          EqBalanceGatheringCurrents, EqEntryIn, RuleEntryOut, 
          RuleExitCurrentsIn, RuleExitValues, EqValueAuxiliaryEdges, 
          OutRules, InRules, EqSwitchingByVertex, EqCompCon, Nlhs, 
          CostArgs, Nrhs, costpluscurrents, EqGeneral, CostArgs2};
        ModuleVarsNames = {"VL", "AM", "EVC", "EVTC", "SC", "BG", 
          "EntranceVertices", "InwardVertices", "ExitVertices", 
          "OutwardVertices", "InEdges", "OutEdges", "AuxiliaryGraph", "FG",
           "EL", "BEL", "FVL", "jargs", "uargs", "AllTransitions", 
          "EntryArgs", "EntryDataAssociation", "ExitCosts", "js", "jvars", 
          "us", "uvars", "jts", "jtvars", "SignedCurrents", 
          "SwitchingCosts", "EqPosJs", "EqPosJts", "EqCurrentCompCon", 
          "EqTransitionCompCon", "NoDeadEnds", 
          "EqBalanceSplittingCurrents", "BalanceSplittingCurrents", 
          "NoDeadStarts", "RuleBalanceGatheringCurrents", 
          "BalanceGatheringCurrents", "EqBalanceGatheringCurrents", 
          "EqEntryIn", "RuleEntryOut", "RuleExitCurrentsIn", 
          "RuleExitValues", "EqValueAuxiliaryEdges", "OutRules", "InRules",
           "EqSwitchingByVertex", "EqCompCon", "Nlhs", "CostArgs", 
          "Nrhs", "costpluscurrents", "EqGeneral", "CostArgs2"};
        Join[Data, AssociationThread[ModuleVarsNames, ModuleVars]]
    ];
NumberVectorQ::usage = 
"NumberVectorQ[j] returns True if the vetor j is numeric."
NumberVectorQ[j_] :=
    And @@ (NumberQ /@ j);

GetKirchhoffMatrix::usage = "GetKirchhoffMatrix[d2e] returns the \
entry current vector, Kirchhoff matrix,  (critical congestion) cost \
function, and the variables in the order corresponding to the Kirchhoff matrix."
GetKirchhoffMatrix[Eqs_] :=
    Module[ {Kirchhoff, EqEntryIn = Lookup[Eqs, "EqEntryIn", True], 
      BalanceGatheringCurrents = 
       Lookup[Eqs, "BalanceGatheringCurrents", {}], 
      BalanceSplittingCurrents = 
       Lookup[Eqs, "BalanceSplittingCurrents", {}], 
      RuleExitCurrentsIn = Lookup[Eqs, "RuleExitCurrentsIn", {}], 
      RuleEntryOut = Lookup[Eqs, "RuleEntryOut", {}], BM, KM, vars, 
      CostArgs = Lookup[Eqs, "CostArgs", <||>], 
      jvars = Lookup[Eqs, "jvars", {}], 
      jtvars = Lookup[Eqs, "jtvars", {}], cost, CCost},
        Kirchhoff = Join[EqEntryIn, (# == 0 & /@ (BalanceGatheringCurrents + BalanceSplittingCurrents))];
        Kirchhoff = Kirchhoff /. Join[RuleExitCurrentsIn, RuleEntryOut];
        vars = Select[Values@Join[jvars,jtvars],MemberQ[Variables[Kirchhoff /. Equal -> List],#]&];
        {BM, KM} = CoefficientArrays[Kirchhoff, vars];
        cost = AssociationThread[vars, KeyMap[Join[jvars,jtvars]][CostArgs]/@vars];
        CCost = cost /@ vars /. MapThread[Rule, {vars, #}] &;
        {-BM, KM, CCost, vars}
    ];
    
MFGPreprocessing::usage =
"MFGPreprocessing[Eqs] returns the association Eqs with the preliminary solution \"InitRules\" and corresponding \'reduced\' \"NewSystem\"."
MFGPreprocessing[Eqs_] :=
    Module[ {InitRules, RuleBalanceGatheringCurrents, EqEntryIn, 
      RuleEntryOut, RuleExitCurrentsIn, RuleExitValues, 
      EqValueAuxiliaryEdges, EqSwitchingByVertex, EqCompCon, 
      EqBalanceSplittingCurrents, EqCurrentCompCon, EqTransitionCompCon,
      EqPosJs, EqPosJts, ModuleVarsNames, ModulesVars, NewSystem, 
      Rules, EqGeneral},
        RuleBalanceGatheringCurrents = Lookup[Eqs, "RuleBalanceGatheringCurrents", $Failed];
        EqEntryIn = Lookup[Eqs, "EqEntryIn", $Failed];
        RuleEntryOut = Lookup[Eqs, "RuleEntryOut", $Failed];
        RuleExitCurrentsIn = Lookup[Eqs, "RuleExitCurrentsIn", $Failed];
        RuleExitValues = Lookup[Eqs, "RuleExitValues", $Failed];
        EqCompCon = Lookup[Eqs, "EqCompCon", $Failed];
        EqCurrentCompCon = Lookup[Eqs, "EqCurrentCompCon", $Failed];
        EqTransitionCompCon = Lookup[Eqs, "EqTransitionCompCon", $Failed];
        EqPosJs = Lookup[Eqs, "EqPosJs", $Failed];
        EqPosJts = Lookup[Eqs, "EqPosJts", $Failed];
        EqGeneral = Lookup[Eqs, "EqGeneral", $Failed];
        EqSwitchingByVertex = Lookup[Eqs, "EqSwitchingByVertex", $Failed];
        EqBalanceSplittingCurrents = Lookup[Eqs, "EqBalanceSplittingCurrents", $Failed];
        EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", $Failed];
        (*First rules:entry currents*)
        InitRules = Association[Flatten[ToRules /@ EqEntryIn]];
        (*no exit at the entrances*)
        AssociateTo[InitRules, RuleEntryOut];
        (*no entrance at the exits*)
        AssociateTo[InitRules, RuleExitCurrentsIn];
        (*currents gathered from transition currents*)
        AssociateTo[InitRules, RuleBalanceGatheringCurrents];
        (*value function:exit costs*)
        AssociateTo[InitRules, RuleExitValues];
        EqSwitchingByVertex = And @@ (Simplify /@ ( EqSwitchingByVertex /. InitRules));
        EqBalanceSplittingCurrents = EqBalanceSplittingCurrents /. InitRules;
        EqValueAuxiliaryEdges = EqValueAuxiliaryEdges /. InitRules;
        Rules = First@Solve[EqBalanceSplittingCurrents && EqValueAuxiliaryEdges] // Quiet;
        InitRules = Join[InitRules /. Rules, Association@Rules];
        NewSystem = 
         MapThread[
          And, {{True, EqPosJts && EqPosJs, EqCurrentCompCon && EqTransitionCompCon && EqCompCon}, 
           Sys2Triple[EqSwitchingByVertex]}
         ];
        {NewSystem, InitRules} = TripleClean[{NewSystem, InitRules}];
        NewSystem[[1]] = EqGeneral;
        {NewSystem, InitRules} = TripleClean[{NewSystem, InitRules}];
        (*sometimes Reduce does not reduce, but if we simplify each term it does.*)
        NewSystem[[3]] = Simplify /@ NewSystem[[3]];
        ModuleVarsNames = {"InitRules", "NewSystem"};
        ModulesVars = {InitRules, NewSystem};
        Join[Eqs, AssociationThread[ModuleVarsNames, ModulesVars]]
    ];

CriticalCongestionSolver::usage = 
  "CriticalCongestionSolver[Eqs] returns Eqs with an association \"AssoCritical\" with rules to 
the solution to the critical congestion case";
CriticalCongestionSolver[$Failed] :=
    $Failed

CriticalCongestionSolver[Eqs_] :=
    Module[ {PreEqs, js, AssoCritical, time},
        {time, PreEqs} = AbsoluteTiming@MFGPreprocessing[Eqs];
        Print["It took ", time, " seconds to preprocess."];
        js = Lookup[PreEqs, "js",$Failed];
        AssoCritical = MFGSystemSolver[PreEqs][AssociationThread[js, 0 js]];
        Join[PreEqs, Association["AssoCritical"-> AssoCritical]]
    ];

(*CriticalCongestionReduce::usage = 
  "CriticalCongestionReduce[Eqs] is an attempt to a worst solver than CriticalCongestionSolver is.";
CriticalCongestionReduce[Eqs_] :=
    Module[ {PreEqs, js, AssoCritical, time},
        {time, PreEqs} = AbsoluteTiming@MFGPreprocessing[Eqs];
        Print["It took ", time, " seconds to preprocess."];
        js = Lookup[PreEqs, "js",$Failed];
        AssoCritical = MFGSystemReduce[PreEqs][AssociationThread[js, 0 js]];
        Join[PreEqs, Association["AssoCritical"-> AssoCritical]]
    ];*)

MFGSystemSolver::usage = 
"MFGSystemSolver[Eqs][edgeEquations] returns the
association with rules to the solution";
MFGSystemSolver[Eqs_][approxJs_] :=
    Module[ {NewSystem, InitRules, pickOne, vars, System, Ncpc,
        costpluscurrents, us, js, jts, usR, jjtsR},
        us = Lookup[Eqs, "us", $Failed];
        js = Lookup[Eqs, "js", $Failed];
        jts = Lookup[Eqs, "jts", $Failed];
        InitRules = Lookup[Eqs, "InitRules", $Failed];
        NewSystem = Lookup[Eqs, "NewSystem", $Failed];
        costpluscurrents = Lookup[Eqs, "costpluscurrents", $Failed];
        Ncpc = RoundValues @ (Expand/@(costpluscurrents /. approxJs));
        InitRules = Expand/@(InitRules /. Ncpc);
        NewSystem = NewSystem /. Ncpc;
        {NewSystem, InitRules} = FinalClean[{NewSystem, InitRules}];
        NewSystem = Sys2Triple[ And @@(Reduce[#,Reals]&/@NewSystem)];
        {NewSystem, InitRules} = FinalClean[{NewSystem, InitRules}];
        System = And @@ NewSystem;
        Which[System === False, 
            Print["MFGSS: There is no solution"], 
            NewSystem[[1]]&&NewSystem[[3]] === True && NewSystem[[2]] =!= True,
            Print["MFGSS: Multiple solutions: ", NewSystem[[2]]//N];
            usR = Select[us, Not[FreeQ[NewSystem, #]] &];
            jjtsR = Select[Join[js, jts], Not[FreeQ[NewSystem[[2]], #]] &];
            vars = Join[usR, jjtsR];
            (*Have to pick one so that all the currents have numerical values*)
            pickOne = Association @ First @ FindInstance[System && And @@ ((# > 0 )& /@ jjtsR), vars, Reals];
            InitRules = Expand /@ Join[InitRules /. pickOne, pickOne];
            Print["\tPicked one value for the variable(s) ", vars, " ", InitRules/@vars//N, " (respectively)"],
         (*System =!=*) True, (*TODO why would we fall here?*)
         Print["System is ", System](*Reducing in the reals should be a good way of simplifying things here.*)
         ];
         (*Print[InitRules];*)
        InitRules = Join[KeyTake[InitRules, us], KeyTake[InitRules, js], KeyTake[InitRules, jts]];
        InitRules
    ];

(*MFGSystemReduce[Eqs_][approxJs_] :=
    Module[ {NewSystem, InitRules, pickOne, vars, System, Ncpc,
        costpluscurrents, us, js, jts, usR, jjtsR},
        us = Lookup[Eqs, "us", $Failed];
        js = Lookup[Eqs, "js", $Failed];
        jts = Lookup[Eqs, "jts", $Failed];
        InitRules = Lookup[Eqs, "InitRules", $Failed];
        NewSystem = Lookup[Eqs, "NewSystem", $Failed];
        costpluscurrents = Lookup[Eqs, "costpluscurrents", $Failed];
        Ncpc = RoundValues @ (Expand/@(costpluscurrents /. approxJs));
        InitRules = Expand/@(InitRules /. Ncpc);
        NewSystem = NewSystem /. Ncpc;
        {NewSystem, InitRules} = FinalReduce[{NewSystem, InitRules}];
        System = And @@ NewSystem;
        Which[System === False, Print["MFGSS: There is no solution"], 
         System =!= True, 
         NewSystem = Reduce[System, Reals];
         {NewSystem, InitRules} = FinalReduce[{Sys2Triple[NewSystem], InitRules}];
         NewSystem = Reduce[And @@ NewSystem, Reals];
         (*not checking if NewSystem is not True...*)
         Print["MFGSS: Multiple solutions: ", NewSystem (*{NewSystem, InitRules}*)];
         usR = Select[us, Not[FreeQ[NewSystem, #]] &];
         jjtsR = Select[Join[js, jts], Not[FreeQ[NewSystem, #]] &];
         vars = Join[usR, jjtsR];
         (*Have to pick one so that all the currents have numerical values*)
         pickOne = Association @ First @
            FindInstance[NewSystem && And @@ ((# > 0 )& /@ jjtsR), vars, Reals];
         InitRules = Expand /@ Join[InitRules /. pickOne, pickOne];
         Print["\tPicked one value for the variable(s) ", vars, " ", InitRules/@vars, " (respectively)"]
        ];
        InitRules = Join[KeyTake[InitRules, us], KeyTake[InitRules, js], KeyTake[InitRules, jts]];
        InitRules
    ];*)

FinalStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, rules_}] :=
    {{EE, NN, OR}, rules}

FinalStep[{{EE_, NN_, OO_}, rules_}] :=
    Module[ {NewSystem, newrules, sorted = True, time, temp},
        {NewSystem, newrules} = TripleClean[{{EE, NN, OO}, rules}];
        If[ NewSystem[[1]] && NewSystem[[3]] === True,
            Return[{NewSystem, newrules}, Module], 
            sorted = RemoveDuplicates[NewSystem[[3]]]
        ];
        temp = PrintTemporary["Iterative DNF convertion..."];
        {time, NewSystem} = 
        AbsoluteTiming[ZAnd[And @@ Take[NewSystem, {1, 2}], sorted]];
        NotebookDelete[temp];
        Print["Iterative DNF convertion took ", time, " seconds to terminate."];
        temp = PrintTemporary["Reducing... ", NewSystem];
        {time,NewSystem} = AbsoluteTiming@Reduce[NewSystem, Reals];
        NotebookDelete[temp];
        Print["Reducing took ", time, " seconds to terminate."];
        NewSystem = Sys2Triple[NewSystem];
        {NewSystem, newrules}
    ];

(*FinalReduceStep[{{EE_, NN_, OO_}, rules_}] :=
    Module[ {NewSystem, newrules, sorted, time, oo, fst},
        {NewSystem, newrules} = TripleClean[{{EE, NN, OO}, rules}];
        oo = Part[NewSystem, 3];
        If[ oo === True,
            sorted = True,
            sorted = SortOp[Part[NewSystem, 3]]
        ];
        fst = And @@ Take[NewSystem, {1, 2}];
        Print["Reducing ", fst && If[ oo === True,
                                      True,
                                      sorted
                                  ]," ..."];
        {time, NewSystem} = AbsoluteTiming[Reduce[fst 
        && If[ Part[NewSystem, 3] === True,
               True,
               sorted
           ], Reals]];
        Print["Reduce took ", time, " seconds to terminate"];
        If[ Head[NewSystem] === Or,
            NewSystem = Sort /@ NewSystem
        ];
        NewSystem = Sys2Triple[NewSystem];
        {NewSystem, newrules}
    ];
   
FinalReduce[{{EE_, NN_, OR_}, rules_}] :=
    (FixedPoint[FinalReduceStep, {{EE, NN, OR}, rules}])*)

FinalClean::usage =
"FinalClean[{{EE,NN,OR},rules}] performs FinalStep, Reduce, and TripleClean on {{EE,NN,OR},rules}.
The result should always be {{True, True, NewOr}, NewRules}"
FinalClean[{{EE_, NN_, OR_}, rules_}] := With[
	{NewSystemRules = FinalStep[{{EE, NN, OR}, rules}]},
	TripleClean[{(Reduce[#, Reals]&/@NewSystemRules[[1]]),NewSystemRules[[2]]}]
];

Sys2Triple::usage =
"Sys2Triple[sys] retrns a triple with equalities, inequalites, and alternatives, respectively."
Sys2Triple[True] = Table[True, 3]

Sys2Triple[False] = Table[False, 3]

Sys2Triple[system_] :=
    Which[
     Head[system] === And, 
      Module[ {groups, EE, OR, NN},
          groups = GroupBy[List @@ system, Head[#] === Equal &];
          EE = And @@ Lookup[groups, True, {}];
          groups = GroupBy[Lookup[groups, False, {}], Head[#] === Or &];
          OR = And @@ Lookup[groups, True, True];
          NN = And @@ Lookup[groups, False, True];
          {EE, NN, OR}
      ], 
     Head[system] === Equal, 
      {system, True, True}, 
     Head[system] === Or, 
         {True, True, system}, 
     True, 
         {True, system, True}];

TripleStep::usage =
"TripleStep[{{EE,NN,OR},Rules}] returns {{True, NewNN, NewOR}, NewRules}, where NewRules contain the solutions to the all equalities found in the system 
after replacing Rules in {EE,NN,OR}."

TripleStep[{{EEs_, NNs_, ORs_}, rules_List}] :=
    TripleStep[{{EEs, NNs, ORs}, Association@rules}]

TripleStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, rules_}] :=
    {{EE, NN, OR}, rules}

TripleStep[{{EE_, NN_?TrueQ, OR_?TrueQ}, rules_Association}] :=
    Module[ {newrules = {}},
        newrules = First@Solve[EE /. rules] // Quiet;
        newrules = Join[rules /. newrules, Association@newrules];
        {{True, NN, OR}, Expand /@ newrules}
    ];

TripleStep[{{EEs_, NNs_, ORs_}, rules_Association}] :=
    Module[ {EE = EEs /. rules, NN = Simplify /@ (NNs /. rules), 
      OR = Simplify /@ (ORs /. rules), NNE, NNO, ORE, ORN,
      newrules = {}},
        NN = Simplify[NN];
        {NNE, NN, NNO} = Sys2Triple[NN];
        {ORE, ORN, OR} = Sys2Triple[OR];
        EE = EE && NNE && ORE;
        If[ EE =!= True && EE =!= False,
            newrules = First@Solve[EE] // Quiet
        ];
        newrules = Expand /@ Join[rules /. newrules, Association@newrules];
        {{EE, NN, OR}, newrules}
    ];

TripleClean::usage =
"TripleClean[{{EE,NN,OR},Rules}] composes TripleStep until it reaches a fixed point, that is, {{True,NewNN,NewOR},NewRules} such that 
replacement of NewRules in NewNN and NewOR do not produce equalities."
TripleClean[{{EE_, NN_, OR_}, rules_}] := FixedPoint[TripleStep, {{EE, NN, OR}, rules}];

ZAnd::usage =
"
ZAnd[xp,xps] returns a system which is equivalent to xp&&xps in disjunctive normal form. 
ZAnd[xp, And[fst, scd] 
ZAnd[xp, eq] returns xp with the solution of eq replaced in it together with eq.
ZAnd[xp, Or[fst,scd]] returns ZAnd[xp, fst]||ZAnd[xp, scd]";

ZAnd[_, False] :=
    False

ZAnd[False, _] :=
    False

ZAnd[xp_, True] :=
    xp

ZAnd[xp_, eq_Equal] := ReZAnd[xp, True, eq]

ZAnd[xp_, And[fst_,rst_]] :=
    If[ Head[fst] === Or,
        RemoveDuplicates@(ReZAnd[Simplify@xp, rst] /@ fst),
        ReZAnd[Simplify@xp, rst, fst]
    ]

ZAnd[xp_, Or[fst_,scd_]] :=
    RemoveDuplicates@(Or@@(ZAnd[Simplify@xp, #] & /@ {fst,scd}))

ZAnd[xp_, leq_] :=
    With[ {ff = Simplify[xp && leq]},
        Print["ZAnd: ", xp&& leq];
        Print["ZAnd: ", ff];
        ff
    ]

(*Operator form of ReZAnd*)
ReZAnd[xp_, rst_] :=
    ReZAnd[xp, rst, #] &

ReZAnd[xp_, rst_, fst_Equal] :=
    With[ {sfst = Simplify[fst]},
        If[ sfst === False,
            False,
            With[ {fsol = First@Solve@sfst},
                ZAnd[Simplify[(xp /. fsol) && fst], ReplaceSolution[rst, fsol]]
            ]
        ]
    ]
  
ReZAnd[xp_, rst_, fst_] :=
    ZAnd[xp && fst, rst]

ReplaceSolution[rst_?BooleanQ, sol_] :=
    rst

ReplaceSolution[rst_, sol_] :=
    With[ {newrst = rst /. sol},
        If[ Head[newrst] === And,
            And[Simplify@First@newrst, Rest@newrst],
            Simplify[newrst]
        ]
    ]

SortOp = ReverseSortBy[Simplify`SimplifyCount]

RemoveDuplicates::usage =
"RemoveDuplicates[xp] sorts and then DeleteDuplicates. 
We need to sort because DeleteDuplicates only deletes identical expressions.
For example (A&&B)||(B&&A) reduce to (A&&B) only after sorting.
"
RemoveDuplicates[xp_And] :=
    DeleteDuplicates[SortOp[xp]];

RemoveDuplicates[xp_Or] :=
    DeleteDuplicates[SortOp[xp]];

RemoveDuplicates[xp_] :=
    xp