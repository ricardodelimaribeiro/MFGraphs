(* Wolfram Language package *)

D2E::usage = 
"D2E[<|\"Vertices List\" -> {1, 2, 3}, \"Adjacency Matrix\" -> {{0, 0, 0}, {1, 0, 1}, {0, 0, 0}}, 
 \"Entrance Vertices and Currents\" -> {{2, I1}}, 
 \"Exit Vertices and Terminal Costs\" -> {{1, U1}, {3, U2}}, 
 \"Switching Costs\" -> {{1, 2, 3, S1}, {3, 2, 1, S2}}|> returns the equations for the stationary mean-field game on the network.]"

Begin["`Private`"]
D2E[Data_Association] :=
    Module[ {BG, EntranceVertices, InwardVertices, InEdges, ExitVertices, 
    OutwardVertices, OutEdges, AuxiliaryGraph, FG, VL, EL, BEL, jargs, 
    js, jvars, FVL, AllTransitions, jts, jtvars, uargs, us, uvars, 
    EqPosCon, EqCurrentCompCon, EqTransitionCompCon, NoDeadEnds, 
    NoDeadStarts, EntryArgs, RulesCriticalCase, EqSwitchingByVertex,
    EntryDataAssociation, ExitCosts, SwitchingCosts, OutRules, InRules, 
    EqSwitchingConditions, EqCompCon, EqValueAuxiliaryEdges, AllOr, 
    SignedCurrents, Nrhs, Nlhs, MinimalTimeRhs, EqAllAll, EqCriticalCase, 
    EqMinimalTime, EqCcs, EqBalanceSplittingCurrents, AllIneq, InitRules, 
    RuleBalanceGatheringCurrents, RuleEntryIn, RuleExitValues, 
    SC = Lookup[Data, "Switching Costs", {}], EqPosJs, TrueEq,
    EqBalanceGatheringCurrents, costpluscurrents, EqNonCritical, RuleNonCritical,
    RuleNonCritical1, RulesCriticalCase1, EqNonCritical1
    },
    (*Checking consistency on the swithing costs*)
        If[ SC =!= {},
            ConsistentSwithingCosts[{a_, b_, c_, S_}] :=
                Module[ {sc = SC, or, de, bounds},
                    or = Cases[sc, {a, b, _, _}];
                    de = Cases[sc, {_, b, c, _}];
                    bounds = Outer[Plus, Last /@ or, Last /@ de] // Flatten;
                    (*And @@ (NonNegative[#-S] &) /@ bounds*)
                    And @@ (#>=S &) /@ bounds
                ];
            EqCcs = ConsistentSwithingCosts /@ SC;
            EqCcs = Simplify[ConsistentSwithingCosts /@ SC];
            (*Print[EqCcs];*)
            If[ AnyTrue[EqCcs,#===False&],
                Print["DataToEquations: Triangle inequalities for switching costs: ", EqCcs];
                Print["DataToEquations: The switching costs are ",Style["incompatible", Red],". \nStopping!"];
                Return[]
            ]
        ];
        
        (*Begin Internal functions for DataToEquations: *)
        IncomingEdges::usage =
            "IncomingEdges[k] all edges \"oriented\" towards the vertex k";
        IncomingEdges[k_] :=
            {k, #1} & /@ IncidenceList[FG,k];
        OutgoingEdges::usage =
            "OutgoingEdges[k] all edges \"oriented\" away from the vertex k";
        OutgoingEdges[k_] :=
            OtherWay /@ ({k, #} & /@ IncidenceList[FG, k]);
        CurrentCompCon::usage =
            "CurrentCompCon[DirectedEdge[a,b]] returns the system of alternatives for currents";
        CurrentCompCon[a_ \[DirectedEdge] b_] :=
            jvars[{a, a \[DirectedEdge] b}] == 0 || jvars[{b, a \[DirectedEdge] b}] == 0;
        CurrentSplitting::usage = 
            "CurrentSplitting[{c,DirectedEdge[a,b]] returns the equations for splitting currents";
        CurrentSplitting[{c_, a_ \[DirectedEdge] b_}] :=
            Select[AllTransitions, (Take[#, 2] == {c, a \[DirectedEdge] b}) &];
        CurrentGathering[{c_, a_ \[DirectedEdge] b_}] :=
            Select[AllTransitions, (Part[#, {1, 3}] == OtherWay[{c, a \[DirectedEdge] b}]) &];
        TransitionCompCon[{v_, edge1_, edge2_}] :=
            jtvars[{v, edge1, edge2}] == 0 || jtvars[{v, edge2, edge1}] == 0;
        ExitValues[a_ \[DirectedEdge] b_] :=
            Total[uvars /@ {{b, DirectedEdge[a, b]}}] == ExitCosts[b];
        ExitRules[a_ \[DirectedEdge] b_] :=
            Total[uvars /@ {{b, DirectedEdge[a, b]}}] -> ExitCosts[b];
        (*Transu: u1<= u2+S(1,2) for agents going from edge 1 to 2 at some vertex*)
        Transu[{v_, edge1_, edge2_}] :=
            uvars[{v,edge1}] <= uvars[{v,edge2}] + SwitchingCosts[{v,edge1,edge2}];
        (*TransuNoSwitch[{v_, edge1_, edge2_}] :=
           uvars[{v,edge2}] -> uvars[{v,edge1}];*)
        Compu[{v_, edge1_, edge2_}] :=
            (jtvars[{v, edge1, edge2}] == 0) || uvars[{v,edge2}]-uvars[{v,edge1}] + SwitchingCosts[{v,edge1,edge2}] == 0;
        (*Default in a, function, corresponds to the classic critical congestion case*)
        a (*[SignedCurrents[edge], edge]:*) =
        Lookup[Data, "a" , Function[{j, edge}, j]];
        (*End*)
        
        (****Graph stuff****)
        BG = AdjacencyGraph[Data["Vertices List"], Data["Adjacency Matrix"], VertexLabels -> "Name", DirectedEdges -> True];
        EntranceVertices = First /@ Data["Entrance Vertices and Currents"];
        ExitVertices = First /@ Data["Exit Vertices and Terminal Costs"];
        Clear["en*"];
        (*InwardVertices defines auxiliary vertices for the entrance vertices*)
        InwardVertices = AssociationThread[EntranceVertices, Symbol["en" <> ToString[#]] & /@ EntranceVertices];
        Clear["ex*"];
        OutwardVertices = AssociationThread[ExitVertices, Symbol["ex" <> ToString[#]] & /@ ExitVertices];
        
        (*InEdges defines auxiliary arguments for the entrance vertices*)
        InEdges = MapThread[DirectedEdge, {InwardVertices /@ EntranceVertices, EntranceVertices}];
        OutEdges = MapThread[DirectedEdge, {ExitVertices, OutwardVertices /@ ExitVertices}];
        AuxiliaryGraph = Graph[Join[InEdges, OutEdges], VertexLabels -> "Name", GraphLayout -> "SpringEmbedding"];
        FG = EdgeAdd[BG, Join[InEdges, OutEdges]];
        VL = Data["Vertices List"];
        EL = EdgeList[FG];
        BEL = EdgeList[BG];
        FVL = VertexList[FG];
        
        
        (*arguments*)
        jargs = Flatten[#, 1] & @ ({AtTail@#, AtHead@#} & /@ EL);
        uargs = jargs;
        AllTransitions = TransitionsAt[FG, #] & /@ FVL // Catenate(*at vertex from first edge to second edge*);
        EntryArgs = AtHead /@ ((EdgeList[AuxiliaryGraph, _ \[DirectedEdge] #] & /@ (First /@ Data["Entrance Vertices and Currents"])) // Flatten[#, 1] &);
        EntryDataAssociation = RoundValues @ AssociationThread[EntryArgs, Last /@ Data["Entrance Vertices and Currents"]];
        ExitCosts = AssociationThread[OutwardVertices /@ (First /@ Data["Exit Vertices and Terminal Costs"]), Last /@ Data["Exit Vertices and Terminal Costs"]];


        (*variables*)
        js = Table[Symbol["j" <> ToString[k]], {k, 1, Length @ jargs}];
        jvars = AssociationThread[jargs, js];
        jts = Table[Symbol["jt" <> ToString[k]], {k, 1, Length @ AllTransitions}];
        jtvars = AssociationThread[AllTransitions, jts];
        us = Table[Symbol["u" <> ToString[k]], {k, 1, Length @ uargs}];
        uvars = AssociationThread[uargs, us];
        costpluscurrents = Table[Symbol["cpc" <> ToString[k]], {k, 1, Length @ BEL}];
        (*costpluscurrentsvars = AssociationThread[BEL, costpluscurrents];*)
        SignedCurrents =   AssociationThread[BEL, (jvars[AtHead[#]] - jvars[AtTail[#]] &) /@ BEL];
        Print["Variables are all set"];
        
        (*Elements of the system*)
        (*Swithing cost is initialized with 0. 
        AssociationThread associates the last association!*)
        SwitchingCosts = AssociationThread[Join[AllTransitions, triple2path[Take[#, 3], FG] & /@ SC], 
            Join[0&/@ AllTransitions, Last[#] & /@ SC]];
        (*Print[SwitchingCosts];*)
        EqPosJs = And @@ ( #>=0& /@ Join[jvars]);(*Inequality*)
        EqPosCon = And @@ ( #>=0& /@ Join[jvars, jtvars]);(*Inequality*)
        EqCurrentCompCon = And @@ (CurrentCompCon /@ EL);(*Or*)
        EqTransitionCompCon = And @@ ((Sort /@ TransitionCompCon /@ AllTransitions) // Union);(*Or*)
        
        NoDeadEnds = IncomingEdges /@ VL // Flatten[#, 1] &;
        EqBalanceSplittingCurrents = And @@ ((jvars[#] == Total[jtvars /@ CurrentSplitting[#]]) & /@ NoDeadEnds);(*Equal*)
        NoDeadStarts = OutgoingEdges /@ VL // Flatten[#, 1] &;
        EqBalanceGatheringCurrents = And @@ ((jvars[#] == Total[jtvars /@ CurrentGathering[#]]) & /@ NoDeadStarts);(*Equal*)
        RuleBalanceGatheringCurrents = (jvars[#] -> Total[jtvars /@ CurrentGathering[#]]) & /@ NoDeadStarts;(*Rule*)
        (*Print[EqBalanceSplittingCurrents&&EqBalanceGatheringCurrents];*)
        
        (*First rules: these have some j in terms of jts*)
        InitRules = Association[RuleBalanceGatheringCurrents];
        
        (*EqEntryIn = And @@ ((jvars[#] == EntryDataAssociation[#]) & /@ (AtHead /@ InEdges));(*Equal*)*)
        RuleEntryIn = (jvars[#] -> EntryDataAssociation[#]) & /@ (AtHead /@ InEdges);(*Rule*)
        (* Not necessary!! InitRules = InitRules /. RuleEntryIn;*)
        AssociateTo[InitRules, RuleEntryIn];
        (*EqExitValues = And @@ (ExitValues /@ IncidenceList[AuxiliaryGraph, OutwardVertices /@ ExitVertices]);(*Equal*)*)
        RuleExitValues = ExitRules /@ IncidenceList[AuxiliaryGraph, OutwardVertices /@ ExitVertices];(*Rule*)
        (*Not necessary!! InitRules = InitRules/.RuleExitValues;*)
        AssociateTo[InitRules, RuleExitValues];
        
        (*Print[Transu/@TransitionsAt[FG, #] & /@ FVL];*)
        
        (*Infinite switching costs here prevent the network from sucking agents from the exits.*)
        OutRules = Rule[#, Infinity] & /@ (Outer[Flatten[{AtTail[#1], #2}] &, OutEdges, EL] // Flatten[#, 1] &);
        InRules = Rule[#, Infinity] & /@ (Outer[{#2[[2]], #1, #2} &, IncidenceList[FG, #] & /@ EntranceVertices, InEdges] // Flatten[#, 2] &);
        AssociateTo[SwitchingCosts, Association[OutRules]];
        (*Print[SwitchingCosts];*)
        AssociateTo[SwitchingCosts, Association[InRules]];
        (*Print[SwitchingCosts];*)
        
        Print["Assembled most elements of the system"];
        (*Switching condition equations*)
        EqSwitchingByVertex = Transu/@TransitionsAt[FG, #] & /@ VL;
        
        
        (*Print[EqSwitchingByVertex];
        Print[Transu/@TransitionsAt[FG, #] & /@ FVL];*)
        
        EqSwitchingByVertex = BooleanConvert[Reduce[Reduce@#, Reals], "CNF"] & /@ EqSwitchingByVertex;
        EqSwitchingByVertex = DeleteCases[EqSwitchingByVertex, True];
        Print["CleanEqualities for the switching conditions on each vertex"];
        {EqSwitchingConditions, InitRules} = CleanEqualities[{EqSwitchingByVertex, InitRules}];
        
        
        (*Print[EqSwitchingConditions];*)
        
        
        EqCompCon = And @@ Compu /@ AllTransitions;(*Or*)
        
        
        (*Print[EqCompCon];*)
        
        
        (*EqCompCon = EqCompCon/.InitRules;*)
        
        
        Print["CleanEqualities for the complementary conditions (given the rules)"];
        {EqCompCon, InitRules} = CleanEqualities[{EqCompCon, InitRules}];
        
        (*Print[EqCompCon];*)
        
        
        
        EqValueAuxiliaryEdges = And @@ ((uvars[AtHead[#]] == uvars[AtTail[#]]) & /@ EdgeList[AuxiliaryGraph]);(*Equal*)
        Print["CleanEqualities for the values at the auxiliary edges"];
        {TrueEq, InitRules} = CleanEqualities[{EqValueAuxiliaryEdges, InitRules}];
        MinimalTimeRhs = Flatten[-a[SignedCurrents[#], #] + SignedCurrents[#] & /@ BEL];
        Nlhs = Flatten[uvars[AtHead[#]] - uvars[AtTail[#]] + SignedCurrents[#] & /@ BEL];
        EqMinimalTime = And @@ (MapThread[(#1 == #2) &, {Nlhs, MinimalTimeRhs}]);
        Nrhs =  Flatten[Cost[SignedCurrents[#], #] + SignedCurrents[#] & /@ BEL];(*one possible cost is IntM*)
        Print["CleanEqualities for the balance conditions in terms of (mostly) transition currents"];
        {TrueEq, InitRules} = CleanEqualities[{EqBalanceSplittingCurrents, InitRules}];
        (*Print[InitRules];Print[EqCompCon/.InitRules];*)
        
        AllOr = EqCurrentCompCon && EqTransitionCompCon && EqCompCon;
        AllOr = AllOr/.InitRules;
        (*Print["Simplifying positivity conditions..."];
        EqPosCon = Simplify[EqPosCon /. InitRules];*)
        EqPosCon = EqPosCon /. InitRules;
        (*Print["Simplifying switching conditions..."];
        EqSwitchingConditions = Simplify[EqSwitchingConditions/.InitRules];*)
        EqSwitchingConditions = EqSwitchingConditions/.InitRules;
        (*Print["Simplifying all inequalities"];
        AllIneq = Simplify[EqPosCon && EqSwitchingConditions];*)
        AllIneq = EqPosCon && EqSwitchingConditions;
        EqAllAll = AllOr && AllIneq;
        EqCriticalCase = And @@ ((# == 0) & /@ Nlhs);(*Equal*)
        (*Print[EqCriticalCase];
        Print[EqCriticalCase/.AssociationThread[us, us/.InitRules]];*)
        
        (*EqNonCritical = And @@ (MapThread[(Equal[#1,#2]) &, {Nlhs, costpluscurrents}])*)
        EqNonCritical = And @@ (MapThread[ Equal[#1,#2] &, {Nlhs, costpluscurrents}])/.AssociationThread[us, us/.InitRules];
        RuleNonCritical1 = Solve[EqNonCritical,Reals]//Quiet;
        
        Print["CleanEqualities for the (non) critical case equations"];
        {TrueEq, RuleNonCritical} = CleanEqualities[{EqNonCritical1, InitRules }];
        (*Print[Expand/@RuleNonCritical];*)
        EqCriticalCase = EqNonCritical /. AssociationThread[costpluscurrents, 0&/@costpluscurrents];
        (*Print[EqCriticalCase];*)
        RulesCriticalCase1 = Expand /@ (RuleNonCritical1 /. AssociationThread[costpluscurrents, 0&/@costpluscurrents]);
        RulesCriticalCase = Expand /@ (RuleNonCritical /. AssociationThread[costpluscurrents, 0&/@costpluscurrents]);
        
        (*Print[Solve[EqCriticalCase, js]]; 
        (*RulesCriticalCaseJs = Association@First[Solve[EqCriticalCase/.RuleEntryIn, js]];*)
        Print["CleanEqualities for the critical case equations"];
        {TrueEq, RulesCriticalCase} = CleanEqualities[{EqCriticalCase, InitRules}];
        
        Print["Expanding critical rules..."];
        RulesCriticalCase = Expand/@RulesCriticalCase;
        Print[RulesCriticalCase];*)
        Join[Data,Association[
        (*Graph structure*)
        "BG" -> BG, 
        "InEdges" -> InEdges, 
        "OutEdges" -> OutEdges,
        "FG" -> FG, 
        (*variables*)
        "jvars" -> jvars, 
        "jtvars" -> jtvars, 
        "uvars" -> uvars,
        "costpluscurrents" -> costpluscurrents,
        "jays" -> SignedCurrents, 
        (*equations*)
        (*complementarity*)
        "AllOr" -> AllOr, (*union of all complementarity conditions*)
        "AllIneq" -> AllIneq,
        "EqPosJs"->EqPosJs,
        "EqCurrentCompCon" -> EqCurrentCompCon,
        "EqTransitionCompCon" -> EqTransitionCompCon,
        (*linear equations (and inequalities)*)
        "EqAllAll" -> EqAllAll,
        "BoundaryRules" -> InitRules,
        "InitRules" -> InitRules,
        "RulesCriticalCase" -> RulesCriticalCase,
        "RulesCriticalCase1" -> RulesCriticalCase1,
        (*"RulesCriticalCaseJs" -> RulesCriticalCaseJs,*)
        "RuleEntryIn" -> RuleEntryIn,
        "Nlhs" -> Nlhs,
        "MinimalTimeRhs" -> MinimalTimeRhs,
        "EqCriticalCase" -> EqCriticalCase,
        "EqNonCritical" ->EqNonCritical,
        "EqNonCritical1" -> EqNonCritical1,
        "RuleNonCritical" -> RuleNonCritical,
        "RuleNonCritical1" -> RuleNonCritical1,
        "EqSwitchingConditions" -> And @@ EqSwitchingByVertex,
		"EqValueAuxiliaryEdges" -> EqValueAuxiliaryEdges,
        "EqBalanceSplittingCurrents" -> EqBalanceSplittingCurrents, 
 		"EqBalanceGatheringCurrents" -> EqBalanceGatheringCurrents,
        "Nrhs" -> Nrhs
        ]
        ]
    ]
    
      
        (*"EntranceVertices" -> EntranceVertices,*) 
        (*"InwardVertices" -> InwardVertices, *)
        (*"ExitVertices" -> ExitVertices,*) 
        (*"OutwardVertices" -> OutwardVertices, *)
        (*"AuxiliaryGraph" -> AuxiliaryGraph,*) 
        (*"VL" -> VL, 
        "EL" -> EL, 
        "BEL" -> BEL, 
        "FVL" -> FVL,*) 
        (*"AllTransitions" -> AllTransitions, 
        "NoDeadEnds" -> NoDeadEnds, 
        "NoDeadStarts" -> NoDeadStarts,*) 
        (*"jargs" -> jargs, 
        "js" -> js,*) 
        (*"jts" -> jts,*) 
        (*"uargs" -> uargs, 
        "us" -> us,*) 
        (*"SwitchingCosts" -> SwitchingCosts,*) 
        (*"OutRules" -> OutRules, *)
        (*"InRules" -> InRules,*) 
        (*"EntryArgs" -> EntryArgs,*) 
        (*"EntryDataAssociation" -> EntryDataAssociation,*) 
        (*"ExitCosts" -> ExitCosts,*) 
        (* 
         
        "EqCompCon" -> EqCompCon,*) 
        (*"EqPosCon" -> EqPosCon, 
        "EqEntryIn" -> EqEntryIn, 
        "EqExitValues" -> EqExitValues, 
        
        "EqValueAuxiliaryEdges" -> EqValueAuxiliaryEdges, 
        "AllIneq" -> AllIneq,*) 
                
      
        
End[]