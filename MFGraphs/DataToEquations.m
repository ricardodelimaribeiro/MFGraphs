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
    NoDeadStarts, EqBalanceSplittingCurrents, EqBalanceGatheringCurrents, EntryArgs, 
    EntryDataAssociation, EqEntryIn, ExitCosts, EqExitValues, SwitchingCosts, OutRules, InRules, 
    EqSwitchingConditions, EqCompCon, EqValueAuxiliaryEdges, AllOr, EqAll, SignedCurrents, Nrhs, Nlhs, MinimalTimeRhs,
    EqAllAll, EqCriticalCase, EqMinimalTime, EqCcs, AllEq, AllIneq, InitRules},
    (*Checking consistency on the swithing costs*)
        If[ Lookup[Data, "Switching Costs", {}] =!= {},
            ConsistentSwithingCosts[{a_, b_, c_, S_}] :=
                Module[ {sc = Data["Switching Costs"], or, de, bounds},
                    or = Cases[sc, {a, b, _, _}];
                    de = Cases[sc, {_, b, c, _}];
                    bounds = Outer[Plus, Last /@ or, Last /@ de] // Flatten;
                    And @@ (NonNegative[#-S] &) /@ bounds
                ];
            EqCcs = Simplify[ConsistentSwithingCosts /@ Data["Switching Costs"]];
            If[ Reduce[EqCcs]==False,
                Print["DataToEquations: Triangle inequalities for switching costs: ", EqCcs];
                Print["DataToEquations: The switching costs are ",Style["incompatible", Red],". \nStopping!"];
                Return[]
            ],
            Print["No swithing costs, value functions are equal at the vertices!"]
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
        
        TransuNoSwitch[{v_, edge1_, edge2_}] :=
           uvars[{v,edge1}] -> uvars[{v,edge2}];
        
        Compu[{v_, edge1_, edge2_}] :=
            (jtvars[{v, edge1, edge2}] == 0) || uvars[{v,edge2}]-uvars[{v,edge1}] + SwitchingCosts[{v,edge1,edge2}] == 0;
        (*Default in a, function, corresponds to the classic critical congestion case*)
        a(*[SignedCurrents[edge], edge]:*) =
        Lookup[Data, "a" , Function[{j, edge}, j]];
        (*End*)
        
        (*Graph stuff*)
        BG = AdjacencyGraph[Data["Vertices List"], Data["Adjacency Matrix"], VertexLabels -> "Name", DirectedEdges -> True];
        EntranceVertices = First /@ Data["Entrance Vertices and Currents"];
        ExitVertices = First /@ Data["Exit Vertices and Terminal Costs"];
        AuxiliaryGraph = Graph[Join[InEdges, OutEdges], VertexLabels -> "Name", GraphLayout -> "SpringEmbedding"];
        FG = EdgeAdd[BG, Join[InEdges, OutEdges]];
        VL = Data["Vertices List"];
        EL = EdgeList[FG];
        BEL = EdgeList[BG];
        FVL = VertexList[FG];
        
        
        Clear["en*"];
        (*InwardVertices defines auxiliary vertices for the entrance vertices*)
        InwardVertices = AssociationThread[EntranceVertices, Symbol["en" <> ToString[#]] & /@ EntranceVertices]; 
        Clear["ex*"];
        OutwardVertices = AssociationThread[ExitVertices, Symbol["ex" <> ToString[#]] & /@ ExitVertices];
        
        (*arguments*)
        (*InEdges defines auxiliary arguments for the entrance vertices*)
        InEdges = MapThread[DirectedEdge, {InwardVertices /@ EntranceVertices, EntranceVertices}];
        OutEdges = MapThread[DirectedEdge, {ExitVertices, OutwardVertices /@ ExitVertices}];
        jargs = Flatten[#, 1] &@({AtTail@#, AtHead@#} & /@ EL);(*Join[AtHead /@ EL, AtTail /@ EL];*)
        uargs = jargs;
        AllTransitions = TransitionsAt[FG, #] & /@ FVL // Catenate(*at vertex from first edge to second edge*);
        BasicTransitions = TransitionsAt[BG, #] & /@ VL // Catenate(*at vertex from first edge to second edge*);
        Print[AllTransitions,"\n",BasicTransitions];
        EntryArgs = AtHead /@ ((EdgeList[AuxiliaryGraph, _ \[DirectedEdge] #] & /@ (First /@ Data["Entrance Vertices and Currents"])) // Flatten[#, 1] &);

        EntryDataAssociation = RoundValues @ AssociationThread[EntryArgs, Last /@ Data["Entrance Vertices and Currents"]];
        ExitCosts = AssociationThread[OutwardVertices /@ (First /@ Data["Exit Vertices and Terminal Costs"]), Last /@ Data["Exit Vertices and Terminal Costs"]];


        (*variables*)
        Clear["j*"];
        js = Table[Symbol["j" <> ToString[k]], {k, 1, Length @ jargs}];
        jvars = AssociationThread[jargs, js];
        Clear["jt*"];
        jts = Table[Symbol["jt" <> ToString[k]], {k, 1, Length @ AllTransitions}];
        jtvars = AssociationThread[AllTransitions, jts];
        Clear["u*"];
        us = Table[Symbol["u" <> ToString[k]], {k, 1, Length @ uargs}];
        uvars = AssociationThread[uargs, us];

        SignedCurrents =   AssociationThread[BEL, (jvars[AtHead[#]] - jvars[AtTail[#]] &) /@ BEL];
        
        
        (*Elements of the system*)
		(*Swithing cost is initialized with 0 AssociationThread associates the last association!*)
        SwitchingCosts = AssociationThread[Join[AllTransitions, triple2path[Take[#, 3], FG] & /@ Data["Switching Costs"]], Join[Table[0, Length[AllTransitions]], Last[#] & /@ Data["Switching Costs"]]];

        (*EqPosCon = And @@ (NonNegative /@ Join[jvars, jtvars])*)
        EqPosCon = And @@ ( #>=0& /@ Join[jvars, jtvars]);(*Inequality*)
        EqCurrentCompCon = And @@ (CurrentCompCon /@ EL);(*Or*)
        EqTransitionCompCon = And @@ ((Sort /@ TransitionCompCon /@ AllTransitions) // Union);(*Or*)

        NoDeadEnds = IncomingEdges /@ VL // Flatten[#, 1] &;
        EqBalanceSplittingCurrents = And @@ ((jvars[#] == Total[jtvars /@ CurrentSplitting[#]]) & /@ NoDeadEnds);(*Equal*)
        RuleBalanceSplittingCurrents = ((jvars[#] -> Total[jtvars /@ CurrentSplitting[#]]) & /@ NoDeadEnds);(*Rule*)
        Print["splitting\n",EqBalanceSplittingCurrents];
		Print[RuleBalanceSplittingCurrents];

        NoDeadStarts = OutgoingEdges /@ VL // Flatten[#, 1] &;
        EqBalanceGatheringCurrents = And @@ ((jvars[#] == Total[jtvars /@ CurrentGathering[#]]) & /@ NoDeadStarts);(*Equal*)
        RuleBalanceGatheringCurrents = ((jvars[#] -> Total[jtvars /@ CurrentGathering[#]]) & /@ NoDeadStarts);(*Rule*)
        Print["gathering\n",EqBalanceGatheringCurrents];
        Print[RuleBalanceGatheringCurrents];


        EqEntryIn = And @@ ((jvars[#] == EntryDataAssociation[#]) & /@ (AtHead /@ InEdges));(*Equal*)
        RuleEntryIn = (jvars[#] -> EntryDataAssociation[#]) & /@ (AtHead /@ InEdges);(*Rule*)
        Print["entry\n",EqEntryIn];
        Print[RuleEntryIn];
        
        (*NonZeroEntryCurrents = And @@ (Positive[EntryDataAssociation[#]] & /@ (AtHead /@ InEdges)); (*useful for the general case...*)*)
        
        EqExitValues = And @@ (ExitValues /@ IncidenceList[AuxiliaryGraph, OutwardVertices /@ ExitVertices]);(*Equal*)
        RuleExitValues = ExitRules /@ IncidenceList[AuxiliaryGraph, OutwardVertices /@ ExitVertices];(*Rule*)
        Print["entry\n",EqExitValues];
        Print[RuleExitValues];
        
        
        (*Infinite switching costs here prevent the network from sucking agents from the exits.*)
        OutRules = Rule[#, Infinity] & /@ (Outer[Flatten[{AtTail[#1], #2}] &, OutEdges, EL] // Flatten[#, 1] &);
        InRules = Rule[#, Infinity] & /@ (Outer[{#2[[2]], #1, #2} &, IncidenceList[FG, #] & /@ EntranceVertices, InEdges] // Flatten[#, 2] &);
        AssociateTo[SwitchingCosts, Association[OutRules]];
        AssociateTo[SwitchingCosts, Association[InRules]];
        EqSwitchingConditions = And @@ Transu /@ AllTransitions;
		Print[TransuNoSwitch /@ BasicTransitions];
        Print["switch\n",EqSwitchingConditions];
        EqSwitchingConditions = Reduce[EqSwitchingConditions, Reals];(*Reducing before joining with other equations*)(*Inequalities*)
        Print["reduced\n",EqSwitchingConditions];
        (*TODO change this to u rules *)
        
        EqCompCon = And @@ Compu /@ AllTransitions;(*Or*)
        EqValueAuxiliaryEdges = And @@ ((uvars[AtHead[#]] - uvars[AtTail[#]] == 0) & /@ EdgeList[AuxiliaryGraph]);(*Equal*)
        RuleValueAuxiliaryEdges =  (uvars[AtHead[#]] -> uvars[AtTail[#]]) & /@ EdgeList[AuxiliaryGraph];(*Rule*)(*TODO check if this is the correct order!*)
        Print[RuleValueAuxiliaryEdges];
        MinimalTimeRhs = Flatten[-a[SignedCurrents[#], #] + SignedCurrents[#] & /@ BEL];
        Nlhs = Flatten[uvars[AtHead[#]] - uvars[AtTail[#]] + SignedCurrents[#] & /@ BEL];
        RuleCriticalCase = Flatten[uvars[AtHead[#]] -> uvars[AtTail[#]]  - SignedCurrents[#] & /@ BEL];(*Rule*)
        (*uvars[AtTail[#]] -> uvars[AtHead[#]]  + SignedCurrents[#]*)
        
        Print[RuleCriticalCase];
        EqCriticalCase = And @@ ((# == 0) & /@ Nlhs);(*Equal*)
        Print[EqCriticalCase];
        EqMinimalTime = And @@ (MapThread[(#1 == #2) &, {Nlhs, MinimalTimeRhs}]);
        Nrhs =  Flatten[IntM[SignedCurrents[#], #] + SignedCurrents[#] & /@ BEL];
        (*EqGeneralCase = And @@ (MapThread[(#1 == #2) &, {Nlhs , Nrhs}]);*)
        AllEq = EqBalanceSplittingCurrents && EqBalanceGatheringCurrents && EqEntryIn && EqExitValues && EqValueAuxiliaryEdges;
        AllOr = EqCurrentCompCon && EqTransitionCompCon && EqCompCon;
        AllIneq = EqPosCon && EqSwitchingConditions;
        InitRules = Join[RuleEntryIn,RuleExitValues];
        EqAll = AllEq && AllIneq;
        EqAllAll = (EqAll && AllOr);
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
        "jays" -> SignedCurrents, 
        (*equations*)
        (*complementarity*)
        "AllOr" -> AllOr, (*union of all complementarity conditions*)
        "AllEq" -> AllEq,
        "AllIneq" -> AllIneq,
        (*linear equations (and inequalities)*)
        "EqAllAll" -> EqAllAll,
        "BoundaryRules" -> InitRules,
        "Nlhs" -> Nlhs,
        "MinimalTimeRhs" -> MinimalTimeRhs,
        "EqCriticalCase" -> EqCriticalCase ,
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
        (*"EqCurrentCompCon" -> EqCurrentCompCon, 
        "EqTransitionCompCon" -> EqTransitionCompCon, 
        "EqCompCon" -> EqCompCon,*) 
        (*"EqPosCon" -> EqPosCon, 
        "EqBalanceSplittingCurrents" -> EqBalanceSplittingCurrents, 
        "EqBalanceGatheringCurrents" -> EqBalanceGatheringCurrents, 
        "EqEntryIn" -> EqEntryIn, 
        "EqExitValues" -> EqExitValues, 
        "EqSwitchingConditions" -> EqSwitchingConditions, 
        "EqValueAuxiliaryEdges" -> EqValueAuxiliaryEdges, 
        "EqAll" -> EqAll,*) 
                
      
        
End[]