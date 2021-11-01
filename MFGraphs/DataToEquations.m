(* Wolfram Language package *)

D2E::usage = 
"D2E[<|\"Vertices List\" -> {1, 2, 3}, \"Adjacency Matrix\" -> {{0, 0, 0}, {1, 0, 1}, {0, 0, 0}}, 
 \"Entrance Vertices and Currents\" -> {{2, I1}}, 
 \"Exit Vertices and Terminal Costs\" -> {{1, U1}, {3, U2}}, 
 \"Switching Costs\" -> {{1, 2, 3, S1}, {3, 2, 1, S2}}|> returns the equations for the stationary mean-field game on the network.]"

Begin["`Private`"]
D2E[Data_Association] :=
    Module[ {BG, FG, AuxiliaryGraph, 
    	VL, FVL, EntranceVertices, InwardVertices, ExitVertices, OutwardVertices, 
     	EL, BEL, OutEdges, InEdges, ExitNeighbors,
     	AllTransitions, 
     	EqCcs, SC = Lookup[Data, "Switching Costs", {}],
     	jargs, js, jvars, jts, jtvars, uargs, us, uvars, SignedCurrents, EntryArgs,
     	NoDeadEnds, NoDeadStarts, ExitCosts, EntryDataAssociation, SwitchingCosts,
    	EqPosCon, EqCurrentCompCon, EqTransitionCompCon, EqSwitchingByVertex,
    	EqSwitchingConditions, EqCompCon, EqValueAuxiliaryEdges, 
    	EqBalanceSplittingCurrents, EqBalanceGatheringCurrents,
    	Nrhs, Nlhs, MinimalTimeRhs, 
    	AllOr, EqAllAll, AllIneq,
		EqCriticalCase, EqMinimalTime,  EqNonCritical, EqPosJs, TrueEq, costpluscurrents, 
		RuleBalanceGatheringCurrents, RuleEntryIn, RuleEntryOut, RuleExitValues, RuleExitCurrentsIn, InitRules, OutRules, 
		InRules, RuleNonCritical, RuleNonCritical1, RulesCriticalCase, RulesCriticalCase1(*, EqNonCritical1*)
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
        (*Swithing cost is initialized with 0. AssociationThread associates the last association!*)
        SwitchingCosts = AssociationThread[Join[AllTransitions, triple2path[Take[#, 3], FG] & /@ SC], 
            Join[0&/@ AllTransitions, Last[#] & /@ SC]];
        EqPosJs = And @@ ( #>=0& /@ Join[jvars]);(*Inequality*)
        EqPosCon = And @@ ( #>=0& /@ Join[jvars, jtvars]);(*Inequality*)
        EqCurrentCompCon = And @@ (CurrentCompCon[jvars] /@ EL);(*Or*)
        EqTransitionCompCon = And @@ ((Sort /@ TransitionCompCon[jtvars] /@ AllTransitions) // Union);(*Or*)
        
        (*Balance Splitting Currents in the full graph*)
		NoDeadEnds = IncomingEdges[FG] /@ VL // Flatten[#, 1] &;
        EqBalanceSplittingCurrents = And @@ ((jvars[#] == Total[jtvars /@ CurrentSplitting[AllTransitions][#]]) & /@ NoDeadEnds);(*Equal*)
		
		(*Gathering currents in the inside of the basic graph*)
        NoDeadStarts = OutgoingEdges[BG] /@ VL // Flatten[#, 1] &;
        RuleBalanceGatheringCurrents = (jvars[#] -> Total[jtvars /@ CurrentGathering[AllTransitions][#]]) & /@ NoDeadStarts;(*Rule*)
        (*First rules: these have some j in terms of jts*)
        InitRules = Association[RuleBalanceGatheringCurrents];
        
        (*get equations for the exit currents at the entry vertices*)
        NoDeadStarts = OutgoingEdges[FG] /@ VL // Flatten[#, 1] &;
        EqBalanceGatheringCurrents = And @@ ((jvars[#] == Total[jtvars /@ CurrentGathering[AllTransitions][#]]) & /@ NoDeadStarts);(*Equal*)
        
        (*Incoming currents*)
		RuleEntryIn = (jvars[#] -> EntryDataAssociation[#]) & /@ (AtHead /@ InEdges);(*Rule*)
        
        (*Outgoing currents at entrances*)
        RuleEntryOut= (jvars[#] -> 0) & /@ (AtTail /@ InEdges);(*Rule*)
        RuleEntryIn = Join[RuleEntryIn, RuleEntryOut];
        (* Not necessary to replace RuleEntryIn in InitRules*)
        AssociateTo[InitRules, RuleEntryIn];

        (*Include Gathering currents information in the rules*)
		{TrueEq, InitRules} = CleanEqualities[{EqBalanceGatheringCurrents/.RuleBalanceGatheringCurrents, InitRules}];

        ExitNeighbors = IncidenceList[AuxiliaryGraph, OutwardVertices /@ ExitVertices];

        (*Incoming currents at the exits are zero*)
        RuleExitCurrentsIn = ExitCurrents[jvars]/@ ExitNeighbors;(*Rule*)
        AssociateTo[InitRules, RuleExitCurrentsIn];
        (*Exit values at exit vertices*)
        RuleExitValues = ExitRules[uvars,ExitCosts] /@ ExitNeighbors;(*Rule*)
        (*Not necessary to replace RuleExitValues in InitRules, there are no us up to now.*)
        AssociateTo[InitRules, RuleExitValues];
        
        EqValueAuxiliaryEdges = And @@ ((uvars[AtHead[#]] == uvars[AtTail[#]]) & /@ EdgeList[AuxiliaryGraph]);(*Equal*)
        Print["CleanEqualities for the values at the auxiliary edges"];
        {TrueEq, InitRules} = CleanEqualities[{EqValueAuxiliaryEdges, InitRules}];
        
        (*Infinite switching costs here prevent the network from sucking agents from the exits.*)
        OutRules = Rule[#, Infinity] & /@ (Outer[Flatten[{AtTail[#1], #2}] &, OutEdges, EL] // Flatten[#, 1] &);
        InRules = Rule[#, Infinity] & /@ (Outer[{#2[[2]], #1, #2} &, IncidenceList[FG, #] & /@ EntranceVertices, InEdges] // Flatten[#, 2] &);
        AssociateTo[SwitchingCosts, Association[OutRules]];
        AssociateTo[SwitchingCosts, Association[InRules]];
        
        Print["Assembled most elements of the system"];
        (*Switching condition equations*)
        EqSwitchingByVertex = Transu[uvars,SwitchingCosts]/@TransitionsAt[FG, #] & /@ VL;
       
        (*EqSwitchingByVertex = BooleanConvert[Reduce[Reduce@#, Reals], "CNF"] & /@ EqSwitchingByVertex;*)(*Maybe this is good with nonzeroswitching costs.*)
        EqSwitchingByVertex = BooleanConvert[Reduce[#, Reals], "CNF"] & /@ EqSwitchingByVertex;
        
        EqSwitchingByVertex = DeleteCases[EqSwitchingByVertex, True];
        
        Print["CleanEqualities for the switching conditions on each vertex"];
        {EqSwitchingConditions, InitRules} = CleanEqualities[{EqSwitchingByVertex, InitRules}];
        
        EqCompCon = And @@ Compu[jtvars,uvars,SwitchingCosts] /@ AllTransitions;(*Or*)
        Print["CleanEqualities for the complementary conditions (given the rules)"];
        {EqCompCon, InitRules} = CleanEqualities[{EqCompCon, InitRules}];
        
        (*Default in a, function, corresponds to the classic critical congestion case*)
        a (*[SignedCurrents[edge], edge]:*) =
        Lookup[Data, "a" , Function[{j, edge}, j]];
        
        MinimalTimeRhs = Flatten[-a[SignedCurrents[#], #] + SignedCurrents[#] & /@ BEL];
        Nlhs = Flatten[uvars[AtHead[#]] - uvars[AtTail[#]] + SignedCurrents[#] & /@ BEL];
        EqMinimalTime = And @@ (MapThread[(#1 == #2) &, {Nlhs, MinimalTimeRhs}]);
        
        Nrhs =  Flatten[-Cost[SignedCurrents[#], #] + SignedCurrents[#] & /@ BEL];(*one possible cost is IntM*)
        Print["CleanEqualities for the balance conditions in terms of (mostly) transition currents"];
        {TrueEq, InitRules} = CleanEqualities[{EqBalanceSplittingCurrents, InitRules}];
        (*Print[TrueEq];*)
        AllOr = EqCurrentCompCon && EqTransitionCompCon && EqCompCon;
        AllOr = AllOr/.InitRules;
        AllOr = BooleanConvert[Simplify /@ AllOr, "CNF"];
        (*Print[InitRules];*)
        {AllOr, InitRules} = CleanEqualities[{AllOr, InitRules}];
        (*Print[InitRules];*)
        EqPosCon = EqPosCon /. InitRules;
        EqSwitchingConditions = EqSwitchingConditions/.InitRules;
        
        AllIneq = EqPosCon && EqSwitchingConditions;
        AllIneq = DeleteDuplicates @ AllIneq;
        EqAllAll = AllOr && AllIneq;
        EqCriticalCase = And @@ ((# == 0) & /@ Nlhs);(*Equal*)
        
        EqNonCritical = And @@ (MapThread[ Equal[#1,#2] &, {Nlhs, costpluscurrents/.InitRules}])/.AssociationThread[us, us/.InitRules];
        RuleNonCritical1 = Solve[EqNonCritical, Reals]//Quiet;
        (*The operations below are commutative!
        Print[RuleNonCritical1/. AssociationThread[costpluscurrents, 0&/@costpluscurrents]];
        Print[Solve[EqNonCritical/. AssociationThread[costpluscurrents, 0&/@costpluscurrents],Reals]//Quiet];*)
        
        (*Print[EqCriticalCase/.AssociationThread[us, us/.InitRules]];
        Print[EqNonCritical/. AssociationThread[costpluscurrents, 0&/@costpluscurrents]];
        Print["right?\n", (EqCriticalCase/.AssociationThread[us, us/.InitRules])/.RuleNonCritical1];
        Print["right?\n", (EqCriticalCase /. AssociationThread[costpluscurrents, 0&/@costpluscurrents]) 
        	/.(RuleNonCritical1 /. AssociationThread[costpluscurrents, 0&/@costpluscurrents])];
        Print[EqNonCritical /. RuleNonCritical1];*)
        
        
        
        Print["CleanEqualities for the (non) critical case equations"];
        {TrueEq, RuleNonCritical} = CleanEqualities[{EqNonCritical, InitRules }];
        (*Print[Expand/@RuleNonCritical];*)
        EqCriticalCase = EqNonCritical /. AssociationThread[costpluscurrents, 0&/@costpluscurrents];
        (*Print[EqCriticalCase];*)
        RulesCriticalCase1 = Expand /@ (RuleNonCritical1 /. AssociationThread[costpluscurrents, 0&/@costpluscurrents]);
        RulesCriticalCase = Expand /@ (RuleNonCritical /. AssociationThread[costpluscurrents, 0&/@costpluscurrents]);
        
        (*Here we change the pourpose of costpluscurrents*)
        costpluscurrents = AssociationThread[costpluscurrents, Nrhs];
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
        "EqPosCon" -> EqPosCon, 
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
        (*"EqNonCritical1" -> EqNonCritical1,*)
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