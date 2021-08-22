(* Wolfram Language package *)

D2E::usage = 
"D2E[<|\"Vertices List\" -> {1, 2, 3}, \"Adjacency Matrix\" -> {{0, 0, 0}, {1, 0, 1}, {0, 0, 0}}, 
 \"Entrance Vertices and Currents\" -> {{2, I1}}, 
 \"Exit Vertices and Terminal Costs\" -> {{1, U1}, {3, U2}}, 
 \"Switching Costs\" -> {{1, 2, 3, S1}, {3, 2, 1, S2}}|> returns the equations for the stationary mean-field game on the network.]"

Begin["`Private`"]
D2E[Dat_?AssociationQ] :=
    Module[ {Data = Dat, BG, EntranceVertices, InwardVertices, InEdges, ExitVertices, 
    OutwardVertices, OutEdges, AuxiliaryGraph, FG, VL, EL, BEL, jargs, 
    js, jvars, FVL, AllTransitions, jts, jtvars, uargs, us, uvars, 
    EqPosCon, EqCurrentCompCon, EqTransitionCompCon, NoDeadEnds, 
    NoDeadStarts, EqBalanceSplittingCurrents, EqBalanceGatheringCurrents, EntryArgs, 
    EntryDataAssociation, EqEntryIn, ExitCosts, EqExitValues, SwitchingCosts, OutRules, InRules, 
    EqSwitchingConditions, EqCompCon, EqValueAuxiliaryEdges, AllOr, EqAll, SignedCurrents, Nrhs, Nlhs, MinimalTimeRhs,
    EqAllAll, EqCriticalCase, EqMinimalTime, EqCcs, AllEq, AllIneq, InitRules},
    (*Checking consistency on the swithing costs*)
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
        ];
          (*Begin Internal functions for DataToEquations: *)
        IncomingEdges[k_] :=
            {k, #1} & /@ IncidenceList[FG,k]; (*all edges "oriented" towards the vertex k*)
        OutgoingEdges[k_] :=
            OtherWay /@ ({k, #} & /@ IncidenceList[FG, k]); (*all edges "oriented" away from the vertex k*)
        CurrentCompCon[a_ \[DirectedEdge] b_] :=
            jvars[{a, a \[DirectedEdge] b}] == 0 || jvars[{b, a \[DirectedEdge] b}] == 0;
        (*CurrentCompCon[edge_] := 
        jvars[AtTail[edge]] || jvars[AtHead[edge]];*)
        CurrentSplitting[{c_, a_ \[DirectedEdge] b_}] :=
            Select[AllTransitions, (Take[#, 2] == {c, a \[DirectedEdge] b}) &];
        CurrentGathering[{c_, a_ \[DirectedEdge] b_}] :=
            Select[AllTransitions, (Part[#, {1, 3}] == OtherWay[{c, a \[DirectedEdge] b}]) &];
        TransitionCompCon[{v_, edge1_, edge2_}] :=
            jtvars[{v, edge1, edge2}] == 0 || jtvars[{v, edge2, edge1}] == 0;
        ExitValues[a_ \[DirectedEdge] b_] :=
            Total[uvars /@ {{b, DirectedEdge[a, b]}}] == ExitCosts[b];
        Transu[{v_, edge1_, edge2_}] :=
            uvars[{v,edge2}]-uvars[{v,edge1}] + SwitchingCosts[{v,edge1,edge2}] >= 0;
        (*NonNegative[uvars[{v,edge2}]-uvars[{v,edge1}] + SwitchingCosts[{v,edge1,edge2}]];(*u1<= u2+S(1,2) for agents going from edge 1 to 2 at some vertex*)*)
        Compu[{v_, edge1_, edge2_}] :=
            (jtvars[{v, edge1, edge2}] == 0) || uvars[{v,edge2}]-uvars[{v,edge1}] + SwitchingCosts[{v,edge1,edge2}] == 0;
        a (*[SignedCurrents[edge], edge]:*)=
        Lookup[Data, "a" , Function[{j, edge}, j]];(*Default corresponds to the classic critical congestion case*)
        (*End*)
        BG = AdjacencyGraph[Data["Vertices List"], Data["Adjacency Matrix"], VertexLabels -> "Name", DirectedEdges -> True];
        EntranceVertices = First /@ Data["Entrance Vertices and Currents"];
        Clear["en*"];
        InwardVertices = AssociationThread[EntranceVertices, Symbol["en" <> ToString[#]] & /@ EntranceVertices]; (*Defines auxiliary vertices for the entrance vertices*)
        InEdges = MapThread[DirectedEdge, {InwardVertices /@ EntranceVertices, EntranceVertices}];
        ExitVertices = First /@ Data["Exit Vertices and Terminal Costs"];
        Clear["ex*"];
        OutwardVertices = AssociationThread[ExitVertices, Symbol["ex" <> ToString[#]] & /@ ExitVertices];
        OutEdges = MapThread[DirectedEdge, {ExitVertices, OutwardVertices /@ ExitVertices}];
        AuxiliaryGraph = Graph[Join[InEdges, OutEdges], VertexLabels -> "Name", GraphLayout -> "SpringEmbedding"];
        FG = EdgeAdd[BG, Join[InEdges, OutEdges]];
        VL = Data["Vertices List"];
        EL = EdgeList[FG];
        BEL = EdgeList[BG];
        jargs = Join[AtHead /@ EL, AtTail /@ EL];
        Clear["j*"];
        js = Table[Symbol["j" <> ToString[k]],{k,1, Length@jargs}];
        jvars = AssociationThread[jargs, js];
        FVL = VertexList[FG];
        SignedCurrents =   AssociationThread[BEL, (jvars[AtHead[#]] - jvars[AtTail[#]] &) /@ BEL];
        AllTransitions = TransitionsAt[FG, #] & /@ FVL // Catenate(*at vertex from first edge to second edge*);
        Clear["jt*"];
        jts = Table[Symbol["jt" <> ToString[k]],{k,1, Length@AllTransitions}];
        jtvars = AssociationThread[AllTransitions, jts];
        uargs = Join[AtTail /@ EL, AtHead /@ EL];
        Clear["u*"];
        us = Table[Symbol["u" <> ToString[k]],{k,1, Length@uargs}];
        uvars = AssociationThread[uargs, us];
        (*EqPosCon = And @@ (NonNegative /@ Join[jvars, jtvars])*)
        EqPosCon = And @@ ( #>=0& /@ Join[jvars, jtvars]);(*Inequality*)
        EqCurrentCompCon = And @@ (CurrentCompCon /@ EL);(*Or*)
        EqTransitionCompCon = And @@ ((Sort /@ TransitionCompCon /@ AllTransitions) // Union);(*Or*)
        NoDeadEnds = IncomingEdges /@ VL // Flatten[#, 1] &;
        NoDeadStarts = OutgoingEdges /@ VL // Flatten[#, 1] &;
        EqBalanceSplittingCurrents = And @@ ((jvars[#] == Total[jtvars /@ CurrentSplitting[#]]) & /@ NoDeadEnds);(*Equal*)
        EqBalanceGatheringCurrents = And @@ ((jvars[#] == Total[jtvars /@ CurrentGathering[#]]) & /@ NoDeadStarts);(*Equal*)
        EntryArgs = AtHead /@ ((EdgeList[AuxiliaryGraph, _ \[DirectedEdge] #] & /@ (First /@ Data["Entrance Vertices and Currents"])) // Flatten[#, 1] &);
        EntryDataAssociation = RoundValues@AssociationThread[EntryArgs, Last /@ Data["Entrance Vertices and Currents"]];
        EqEntryIn = And @@ ((jvars[#] == EntryDataAssociation[#]) & /@ (AtHead /@ InEdges));(*Equal*)
        (*NonZeroEntryCurrents = And @@ (Positive[EntryDataAssociation[#]] & /@ (AtHead /@ InEdges)); (*useful for the general case...*)*)
        ExitCosts = AssociationThread[OutwardVertices /@ (First /@ Data["Exit Vertices and Terminal Costs"]), Last /@ Data["Exit Vertices and Terminal Costs"]];
        EqExitValues = And @@ (ExitValues /@ IncidenceList[AuxiliaryGraph, OutwardVertices /@ ExitVertices]);(*Equal*)
        SwitchingCosts = AssociationThread[Join[AllTransitions, triple2path[Take[#, 3], FG] & /@ Data["Switching Costs"]], Join[Table[0, Length[AllTransitions]], Last[#] & /@ Data["Switching Costs"]]];(*Swithing cost is initialized with 0 AssociationThread associates the last association!*)
        (*Infinite switching costs here prevent the network from sucking agents from the exits.*)
        OutRules = Rule[#, Infinity] & /@ (Outer[Flatten[{AtTail[#1], #2}] &, OutEdges, EL] // Flatten[#, 1] &);
        AssociateTo[SwitchingCosts, Association[OutRules]];
        InRules = Rule[#, Infinity] & /@ (Outer[{#2[[2]], #1, #2} &, IncidenceList[FG, #] & /@ EntranceVertices, InEdges] // Flatten[#, 2] &);
        AssociateTo[SwitchingCosts, Association[InRules]];
        EqSwitchingConditions = Reduce[(And @@ Transu /@ AllTransitions), Reals];(*Reducing before joining with other equations*)(*Inequalities*)
        EqCompCon = And @@ Compu /@ AllTransitions;(*Or*)
        EqValueAuxiliaryEdges = And @@ ((uvars[AtHead[#]] - uvars[AtTail[#]] == 0) & /@ EdgeList[AuxiliaryGraph]);(*Equal*)
        MinimalTimeRhs = Flatten[-a[SignedCurrents[#], #] + SignedCurrents[#] & /@ BEL];
        Nlhs = Flatten[uvars[AtHead[#]] - uvars[AtTail[#]] + SignedCurrents[#] & /@ BEL];
        EqCriticalCase = And @@ ((# == 0) & /@ Nlhs);
        EqMinimalTime = And @@ (MapThread[(#1 == #2) &, {Nlhs, MinimalTimeRhs}]);
        Nrhs =  Flatten[IntM[SignedCurrents[#], #] + SignedCurrents[#] & /@ BEL];
        (*EqGeneralCase = And @@ (MapThread[(#1 == #2) &, {Nlhs , Nrhs}]);*)
        AllEq = EqBalanceSplittingCurrents && EqBalanceGatheringCurrents && EqEntryIn && EqExitValues && EqValueAuxiliaryEdges;
        AllOr = EqCurrentCompCon && EqTransitionCompCon && EqCompCon;
        AllIneq = EqPosCon && EqSwitchingConditions;
        InitRules = Join[ToRules[EqEntryIn],ToRules[EqExitValues]];
        EqAll = AllEq && AllIneq;
        EqAllAll = (EqAll && AllOr);
        Join[Data,Association[
        (*Graph structure*)
        "BG" -> BG, 
        (*"EntranceVertices" -> EntranceVertices, 
        "InwardVertices" -> InwardVertices, 
        "InEdges" -> InEdges, 
        "ExitVertices" -> ExitVertices, 
        "OutwardVertices" -> OutwardVertices, 
        "OutEdges" -> OutEdges,*) 
        (*"AuxiliaryGraph" -> AuxiliaryGraph,*) 
        "FG" -> FG, 
        (*"VL" -> VL, 
        "EL" -> EL, 
        "BEL" -> BEL, 
        "FVL" -> FVL,*) 
        (*"AllTransitions" -> AllTransitions, 
        "NoDeadEnds" -> NoDeadEnds, 
        "NoDeadStarts" -> NoDeadStarts,*) 
        (*variables*)
        (*"jargs" -> jargs, 
        "js" -> js,*) 
        "jvars" -> jvars, 
        (*"jts" -> jts,*) 
        "jtvars" -> jtvars, 
        (*"uargs" -> uargs, 
        "us" -> us,*) 
        "uvars" -> uvars,
        (*"SwitchingCosts" -> SwitchingCosts,*) 
        "OutRules" -> OutRules, 
        "InRules" -> InRules, 
        (*"EntryArgs" -> EntryArgs,*) 
        (*"EntryDataAssociation" -> EntryDataAssociation,*) 
        "jays" -> SignedCurrents, 
        "ExitCosts" -> ExitCosts, 
        (*equations*)
        (*complementarity*)
        (*"EqCurrentCompCon" -> EqCurrentCompCon, 
        "EqTransitionCompCon" -> EqTransitionCompCon, 
        "EqCompCon" -> EqCompCon,*) 
        "AllOr" -> AllOr, (*union of all complementarity conditions*)
        "AllEq" -> AllEq,
        "AllIneq" -> AllIneq,
        (*linear equations (and inequalities)*)
        (*"EqPosCon" -> EqPosCon, 
        "EqBalanceSplittingCurrents" -> EqBalanceSplittingCurrents, 
        "EqBalanceGatheringCurrents" -> EqBalanceGatheringCurrents, 
        "EqEntryIn" -> EqEntryIn, 
        "EqExitValues" -> EqExitValues, 
        "EqSwitchingConditions" -> EqSwitchingConditions, 
        "EqValueAuxiliaryEdges" -> EqValueAuxiliaryEdges, 
        "EqAll" -> EqAll,*) 
        "EqAllAll" -> EqAllAll,
        "BoundaryRules" -> InitRules,
        "Nlhs" -> Nlhs,
        "MinimalTimeRhs" -> MinimalTimeRhs,
        "EqCriticalCase" -> EqCriticalCase ,
        "Nrhs" -> Nrhs
        ]
        ]
    ]

End[]