(* Wolfram Language package *)

(*Assembles the stationary mfg equations using the notion of current*)
(*DataToEquations::usage = 
"DataToEquations[<|\"Vertices List\" -> {1, 2, 3}, \"Adjacency Matrix\" -> {{0, 0, 0}, {1, 0, 1}, {0, 0, 0}}, 
 \"Entrance Vertices and Currents\" -> {{2, I1}}, 
 \"Exit Vertices and Terminal Costs\" -> {{1, U1}, {3, U2}}, 
 \"Switching Costs\" -> {{1, 2, 3, S1}, {3, 2, 1, S2}}|>] returns the equations for the stationary mean-field game on the network. 
 This also returns the solution to the critical congestion case and a reduced version of the system/rules for the \"no case\" system.";*)

D2E::usage = 
"D2E[<|\"Vertices List\" -> {1, 2, 3}, \"Adjacency Matrix\" -> {{0, 0, 0}, {1, 0, 1}, {0, 0, 0}}, 
 \"Entrance Vertices and Currents\" -> {{2, I1}}, 
 \"Exit Vertices and Terminal Costs\" -> {{1, U1}, {3, U2}}, 
 \"Switching Costs\" -> {{1, 2, 3, S1}, {3, 2, 1, S2}}|> returns the equations for the stationary mean-field game on the network.]"

(*CriticalBundle::usage =
"CriticalBundle[Data] builds equations and solves the critical congestion case."*)

CriticalCongestionSolver::usage = 
"CriticalCongestionSolver[eq_Association] returns the critical congestion solution."

CriticalCongestionStep::usage =
""
(*CriticalCongestionSolverN::usage = 
""*)

Begin["`Private`"]
D2E[Dat_?AssociationQ] :=
    Module[ {showAll = True, Data = (*RoundValues @ *) Dat, BG, EntranceVertices, InwardVertices, InEdges, ExitVertices, 
        OutwardVertices, OutEdges, AuxiliaryGraph, FG, VL, EL, BEL, jargs, 
        js, jvars, FVL, AllTransitions, jts, jtvars, uargs, us, uvars, 
        EqPosCon, EqCurrentCompCon, EqTransitionCompCon, NoDeadEnds, 
        NoDeadStarts, EqBalanceSplittingCurrents, EqBalanceGatheringCurrents, EntryArgs, 
        EntryDataAssociation, EqEntryIn, NonZeroEntryCurrents, ExitCosts, EqExitValues, SwitchingCosts, OutRules, InRules, 
        EqSwitchingConditions, EqCompCon, EqValueAuxiliaryEdges, EqAllComp, EqAll, SignedCurrents, Nrhs, Nlhs, RulesEntryIn, MinimalTimeRhs,
        RulesExitValues, EqAllRules, EqAllCompRules, EqAllAll, EqAllAllRules, 
        EqCriticalCase, EqMinimalTime, BoundaryRules, EqGeneralCase, EqCcs},
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
        
        (*EqAllComp = EqCurrentCompCon && EqTransitionCompCon && EqCompCon;
        EqAll = EqPosCon && EqBalanceSplittingCurrents && EqBalanceGatheringCurrents && 
            EqEntryIn && EqExitValues && EqValueAuxiliaryEdges && EqSwitchingConditions && 
            And @@ EqCcs;(*includes transition inequalities for the values and the triangle inequalities, too.*)*)
        (*Pre-processing: substitute rules in all equations and hand-sides...*)
        (*EqAllAll = EqAll && EqAllComp;*)
        EqAllAll = (EqAll && AllOr);
        If[ showAll == True,
            Association[{
                "Data" -> Data,
            (*Graph structure*)
            "BG" -> BG, 
            "EntranceVertices" -> EntranceVertices, 
            "InwardVertices" -> InwardVertices, 
            "InEdges" -> InEdges, 
            "ExitVertices" -> ExitVertices, 
            "OutwardVertices" -> OutwardVertices, 
            "OutEdges" -> OutEdges, 
            "AuxiliaryGraph" -> AuxiliaryGraph, 
            "FG" -> FG, 
            "VL" -> VL, 
            "EL" -> EL, 
            "BEL" -> BEL, 
            "FVL" -> FVL, 
            "AllTransitions" -> AllTransitions, 
            "NoDeadEnds" -> NoDeadEnds, 
            "NoDeadStarts" -> NoDeadStarts, 
            (*variables*)
            "jargs" -> jargs, 
            "js" -> js, 
            "jvars" -> jvars, 
            "jts" -> jts, 
            "jtvars" -> jtvars, 
            "uargs" -> uargs, 
            "us" -> us, 
            "uvars" -> uvars,
            "SwitchingCosts" -> SwitchingCosts, 
            "OutRules" -> OutRules, 
            "InRules" -> InRules, 
            
            "EntryArgs" -> EntryArgs, 
            "EntryDataAssociation" -> EntryDataAssociation, 
            "jays" -> SignedCurrents, 
            (*"NonZeroEntryCurrents" -> NonZeroEntryCurrents,*) 
            "ExitCosts" -> ExitCosts, 
            
            (*equations*)
            (*complementarity*)
            "EqCurrentCompCon" -> EqCurrentCompCon, 
            "EqTransitionCompCon" -> EqTransitionCompCon, 
            "EqCompCon" -> EqCompCon, 
            "EqAllComp" -> AllOr, (*union of all complementarity conditions*)
            
            (*linear equations (and inequalities)*)
            "EqPosCon" -> EqPosCon, 
            "EqBalanceSplittingCurrents" -> EqBalanceSplittingCurrents, 
            "EqBalanceGatheringCurrents" -> EqBalanceGatheringCurrents, 
            "EqEntryIn" -> EqEntryIn, 
            "EqExitValues" -> EqExitValues, 
            "EqSwitchingConditions" -> EqSwitchingConditions, 
            "EqValueAuxiliaryEdges" -> EqValueAuxiliaryEdges, 
            "EqAll" -> EqAll, 
            (*"EqAllCompRules" -> EqAllCompRules,*)
            (*"EqAllRules" -> EqAllRules,*)
            "EqAllAll" -> EqAllAll,
            "BoundaryRules" -> InitRules,
            (*"reduced1" -> reduced1,
            "reduced2" -> reduced2,*)
            (*"EqAllAllSimple" -> EqAllAllSimple,*)
            (*"EqAllAllRules" -> EqAllAllRules,*)
            "Nlhs" -> Nlhs,
            "MinimalTimeRhs" -> MinimalTimeRhs,
            "EqCriticalCase" -> EqCriticalCase ,
            (*"criticalreduced2" -> criticalreduced2,*)
            "Nrhs" -> Nrhs(*,
            "EqGeneralCase" -> EqGeneralCase*)
            }],
            Association[{
              "Data" -> Data,
              "FG" -> FG, 
              "EqAllAll" -> EqAllAll,
              "BoundaryRules" -> BoundaryRules,
              "Nlhs" -> Nlhs,
              "EqCriticalCase" -> EqCriticalCase ,
              "Nrhs" -> Nrhs,
              "EqGeneralCase" -> EqGeneralCase
              }]
        ]
    ]



DataToEquations[Data_?AssociationQ] :=
    Module[ {showAll = True, BG, EntranceVertices, InwardVertices, InEdges, ExitVertices, 
    OutwardVertices, OutEdges, AuxiliaryGraph, FG, VL, EL, BEL, jargs, 
    js, jvars, FVL, AllTransitions, jts, jtvars, uargs, us, uvars, 
    EqPosCon, EqCurrentCompCon, EqTransitionCompCon, NoDeadEnds, 
    NoDeadStarts, EqBalanceSplittingCurrents, EqBalanceGatheringCurrents, EntryArgs, 
    EntryDataAssociation, EqEntryIn, NonZeroEntryCurrents, ExitCosts, EqExitValues, SwitchingCosts, OutRules, InRules, 
    EqSwitchingConditions, EqCompCon, EqValueAuxiliaryEdges, EqAllComp, EqAll, SignedCurrents, Nrhs, Nlhs, RulesEntryIn, MinimalTimeRhs,
    RulesExitValues, EqAllRules, EqAllCompRules, EqAllAll, EqAllAllRules, reduced1,reduced2,
    EqCriticalCase, EqMinimalTime, BoundaryRules, criticalreduced1, criticalreduced2, EqGeneralCase, EqCcs, (*EqAllAllSimple,*) sol, system, rules ,TOL, time,nocasesystem,nocaserules(*, result*)},
    
    (*Set the tolerance for Chop*)
        TOL = 10^(-6);
        (*Checking consistency on the swithing costs*)
        ConsistentSwithingCosts[{a_, b_, c_, S_}] :=
            Module[ {sc = Data["Switching Costs"], or, de, bounds},
                or = Cases[sc, {a, b, _, _}];
                de = Cases[sc, {_, b, c, _}];
                bounds = Outer[Plus, Last /@ or, Last /@ de] // Flatten;
                (*And @@ (NonNegative[#-S] &) /@ bounds*)
                And @@ (# >= S &) /@ bounds
            ];
        EqCcs = ConsistentSwithingCosts /@ Data["Switching Costs"];
        Print["DataToEquations: Triangle inequalities for switching costs: ", EqCcs,
            "\nDataToEquations: Reduced: ",Reduce[EqCcs]];
        If[ And @@ EqCcs,
              (*this happens when we have non-false condition.*)
            Print["DataToEquations: The switching costs are compatible."],
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
            NonNegative[uvars[{v,edge2}]-uvars[{v,edge1}] + SwitchingCosts[{v,edge1,edge2}]];(*u1<= u2+S(1,2) for agents going from edge 1 to 2 at some vertex*)
        Compu[{v_, edge1_, edge2_}] :=
            (jtvars[{v, edge1, edge2}] == 0) || uvars[{v,edge2}]-uvars[{v,edge1}] + SwitchingCosts[{v,edge1,edge2}] == 0;
        a (*[SignedCurrents[edge], edge]:*)=
            Lookup[Data, "a" , Function[{j, edge}, j]];(*Default corresponds to the classic critical congestion case*)
        (*End*)
        BG = AdjacencyGraph[Data["Vertices List"], Data["Adjacency Matrix"], VertexLabels -> "Name", DirectedEdges -> True, GraphLayout -> "SpringEmbedding"];
        EntranceVertices = First /@ Data["Entrance Vertices and Currents"];
        InwardVertices = AssociationThread[EntranceVertices, Unique["en"] & /@ EntranceVertices]; (*Defines auxiliary vertices for the entrance vertices*)
        InEdges = MapThread[DirectedEdge, {InwardVertices /@ EntranceVertices, EntranceVertices}];
        ExitVertices = First /@ Data["Exit Vertices and Terminal Costs"];
        OutwardVertices = AssociationThread[ExitVertices, Unique["ex"] & /@ ExitVertices];
        OutEdges = MapThread[DirectedEdge, {ExitVertices, OutwardVertices /@ ExitVertices}];
        AuxiliaryGraph = Graph[Join[InEdges, OutEdges], VertexLabels -> "Name", GraphLayout -> "SpringEmbedding"];
        FG = EdgeAdd[BG, Join[InEdges, OutEdges]];
        VL = Data["Vertices List"];
        EL = EdgeList[FG];
        BEL = EdgeList[BG];
        jargs = Join[AtHead /@ EL, AtTail /@ EL];
        js = Unique["j"] & /@ jargs;
        jvars = AssociationThread[jargs, js];
        FVL = VertexList[FG];
        SignedCurrents =   AssociationThread[BEL, (jvars[AtHead[#]] - jvars[AtTail[#]] &) /@ BEL];
        AllTransitions = TransitionsAt[FG, #] & /@ FVL // Catenate(*at vertex from first edge to second edge*);
        jts = Unique["jt"] & /@ AllTransitions;
        jtvars = AssociationThread[AllTransitions, jts];
        uargs = Join[AtTail /@ EL, AtHead /@ EL];
        us = Unique["u"] & /@ uargs;
        uvars = AssociationThread[uargs, us];
        EqPosCon = And @@ (#>=0& /@ Join[jvars, jtvars]);
        (*EqPosCon = And @@ (NonNegative /@ Join[jvars, jtvars]);*)
        EqCurrentCompCon = And @@ (CurrentCompCon /@ EL);
        EqTransitionCompCon = And @@ ((Sort /@ TransitionCompCon /@ AllTransitions) // Union);
        NoDeadEnds = IncomingEdges /@ VL // Flatten[#, 1] &;
        NoDeadStarts = OutgoingEdges /@ VL // Flatten[#, 1] &;
        EqBalanceSplittingCurrents = And @@ ((jvars[#] == Total[jtvars /@ CurrentSplitting[#]]) & /@ NoDeadEnds);
        EqBalanceGatheringCurrents = And @@ ((jvars[#] == Total[jtvars /@ CurrentGathering[#]]) & /@ NoDeadStarts);
        EqAll = EqPosCon && EqBalanceSplittingCurrents && EqBalanceGatheringCurrents;(*takes too long to Reduce, so don't!*)
        EntryArgs = AtHead /@ ((EdgeList[AuxiliaryGraph, _ \[DirectedEdge] #] & /@ (First /@ Data["Entrance Vertices and Currents"])) // Flatten[#, 1] &);
        EntryDataAssociation = AssociationThread[EntryArgs, Last /@ Data["Entrance Vertices and Currents"]];
        EqEntryIn = And @@ ((jvars[#] == EntryDataAssociation[#]) & /@ (AtHead /@ InEdges));
        RulesEntryIn = ToRules[EqEntryIn];
        NonZeroEntryCurrents = And @@ (Positive[EntryDataAssociation[#]] & /@ (AtHead /@ InEdges)); (*useful for the general case...*)
        EqAll = EqAll && EqEntryIn;
        ExitCosts = AssociationThread[OutwardVertices /@ (First /@ Data["Exit Vertices and Terminal Costs"]), Last /@ Data["Exit Vertices and Terminal Costs"]];
        EqExitValues = And @@ (ExitValues /@ IncidenceList[AuxiliaryGraph, OutwardVertices /@ ExitVertices]);
        RulesExitValues = ToRules[EqExitValues];
        EqAll = EqAll && EqExitValues;
        SwitchingCosts = AssociationThread[Join[AllTransitions, triple2path[Take[#, 3], FG] & /@ Data["Switching Costs"]], Join[Table[0, Length[AllTransitions]], Last[#] & /@ Data["Switching Costs"]]];(*Swithing cost is initialized with 0 AssociationThread associates the last association!*)
        (*Infinite switching costs here prevent the network from sucking agents from the exits.*)
        OutRules = Rule[#, Infinity] & /@ (Outer[Flatten[{AtTail[#1], #2}] &, OutEdges, EL] // Flatten[#, 1] &);
        AssociateTo[SwitchingCosts, Association[OutRules]];
        InRules = Rule[#, Infinity] & /@ (Outer[{#2[[2]], #1, #2} &, IncidenceList[FG, #] & /@ EntranceVertices, InEdges] // Flatten[#, 2] &);
        AssociateTo[SwitchingCosts, Association[InRules]];
        EqSwitchingConditions = Reduce[(And @@ Transu /@ AllTransitions), Reals];(*Reducing before joining with other equations*)
        EqAll = EqAll && EqSwitchingConditions && And @@ EqCcs;(*includes transition inequalities for the values and the triangle inequalities, too.*)
        EqCompCon = And @@ Compu /@ AllTransitions;
        EqValueAuxiliaryEdges = And @@ ((uvars[AtHead[#]] - uvars[AtTail[#]] == 0) & /@ EdgeList[AuxiliaryGraph]);
        EqAllComp = EqCurrentCompCon && EqTransitionCompCon && EqCompCon;
        EqAll = EqAll && EqValueAuxiliaryEdges;
        (*Pre-processing: substitute rules in all equations and hand-sides...*)
        BoundaryRules = Join[RulesEntryIn, RulesExitValues];
        EqAllAll = EqAll && EqAllComp;
        EqAllAllRules = EqAllCompRules && EqAllRules;
        Print["DataToEquations: Finished assembling strucural equations."];
        Print["Reducing the structural system ... "];
        MinimalTimeRhs = Flatten[-a[SignedCurrents[#], #] + SignedCurrents[#] & /@ BEL];
        Nlhs = Flatten[uvars[AtHead[#]] - uvars[AtTail[#]] + SignedCurrents[#] & /@ BEL];
        Print["DataToEquations: Critical case ... "];
        EqCriticalCase = And @@ ((# == 0) & /@ Nlhs);
        EqMinimalTime = And @@ (MapThread[(#1 == #2) &, {Nlhs, MinimalTimeRhs}]);
        {system, rules} = CleanEqualities[{(EqAllAll && EqCriticalCase), {}}];
        {nocasesystem, nocaserules} = CleanEqualities[{EqAllAll, {}}];
        {time,nocasesystem} = AbsoluteTiming @ NewReduce[nocasesystem]; (*does this step help speed up fixedreducex1?*)
        Print["DataToEquations: It took ", time, " seconds to reduce the general case with NewReduce.","\nThe system is :", nocasesystem];
 
        (*Print["DataToEquations: The system is:\n", system,
            "\nand the rules are:\n", rules];*)
            
            (*There could be some alternatives and inequalities left, so we BooleanConvert and solve/replace throught the tree*)
        {time,system} = AbsoluteTiming @ NewReduce[system];
        Print["DataToEquations: It took ", time, " seconds to reduce the critical congestion case with NewReduce!","\nThe system is ", system];
        (*Print["DataToEquations: After NewReduce, the system is:\n", system,
            "\nand the rules are:\n", rules];*)
        Which[system === True ||
            Head[system] === Equal || 
            (Head[system] === And && DeleteDuplicates[Head /@ system] === Equal),
            sol = Solve[system, Reals];
            If[ Length[sol] == 1,
                {system, rules} = {system, AssociateTo[ rules, First@ sol]} /. First@ sol,
                Print[sol];
            ],
            system === False,
                Print["DataToEquations: ",Style["Incompatible system", Red],"."],
            True,
            Print["DataToEquations: The system does not have the original structure: ", system];
            {system,rules} = CleanEqualities[{system, rules}];
            system = Reduce[system, Reals];
            Print["DataToEquations: ", Style["Multiple solutions", Red]," \n\t\t", system, "\n\tFinding one instance..."];
            (*FindInstance!!!*)
            Module[ {solved = And @@ (Outer[Equal,Keys[rules],Values[rules]]//Diagonal) , onesol},
                onesol = First@FindInstance[system && solved, Join[js, us, jts], Reals];    
                (*TODO: should we choose a method to select a solution?*)
                rules = Association[onesol];
                system = system/.rules
            ];
        ];
        criticalreduced1 = {system, Simplify /@ rules};
        If[ system === True,
            Print["DataToEquations: Critical congestion solved."];,
            Print["DataToEquations: There are ", Style["multiple solutions", Red]," for the given data."];
        ];
        Nrhs =  Flatten[IntM[SignedCurrents[#], #] + SignedCurrents[#] & /@ BEL];
        EqGeneralCase = And @@ (MapThread[(#1 == #2) &, {Nlhs , Nrhs}]);
        
        
        (*Print["DataToEquations: Done."];*)
        If[ showAll == True,
            Association[{
                "Data" -> Data,
            (*Graph structure*)
            "BG" -> BG, 
            "EntranceVertices" -> EntranceVertices, 
            "InwardVertices" -> InwardVertices, 
            "InEdges" -> InEdges, 
            "ExitVertices" -> ExitVertices, 
            "OutwardVertices" -> OutwardVertices, 
            "OutEdges" -> OutEdges, 
            "AuxiliaryGraph" -> AuxiliaryGraph, 
            "FG" -> FG, 
            "VL" -> VL, 
            "EL" -> EL, 
            "BEL" -> BEL, 
            "FVL" -> FVL, 
            "AllTransitions" -> AllTransitions, 
            "NoDeadEnds" -> NoDeadEnds, 
            "NoDeadStarts" -> NoDeadStarts, 
            (*variables*)
            "jargs" -> jargs, 
            "js" -> js, 
            "jvars" -> jvars, 
            "jts" -> jts, 
            "jtvars" -> jtvars, 
            "uargs" -> uargs, 
            "us" -> us, 
            "uvars" -> uvars,
            "SwitchingCosts" -> SwitchingCosts, 
            "OutRules" -> OutRules, 
            "InRules" -> InRules, 
            
            
            
            "EntryArgs" -> EntryArgs, 
            "EntryDataAssociation" -> EntryDataAssociation, 
            "jays" -> SignedCurrents, 
            "NonZeroEntryCurrents" -> NonZeroEntryCurrents, 
            "ExitCosts" -> ExitCosts, 
            
            
            
            (*equations*)
            (*complementarity*)
            "EqCurrentCompCon" -> EqCurrentCompCon, 
            "EqTransitionCompCon" -> EqTransitionCompCon, 
            "EqCompCon" -> EqCompCon, 
            "EqAllComp" -> EqAllComp, (*union of all complementarity conditions*)
            
            (*linear equations (and inequalities)*)
            "EqPosCon" -> EqPosCon, 
            "EqBalanceSplittingCurrents" -> EqBalanceSplittingCurrents, 
            "EqBalanceGatheringCurrents" -> EqBalanceGatheringCurrents, 
            "EqEntryIn" -> EqEntryIn, 
            "EqExitValues" -> EqExitValues, 
            "EqSwitchingConditions" -> EqSwitchingConditions, 
            "EqValueAuxiliaryEdges" -> EqValueAuxiliaryEdges, 
            "EqAll" -> EqAll, 
            (*"EqAllCompRules" -> EqAllCompRules,*)
            (*"EqAllRules" -> EqAllRules,*)
            "EqAllAll" -> EqAllAll,
            "BoundaryRules" -> BoundaryRules,
            (*"reduced1" -> reduced1,
            "reduced2" -> reduced2,*)
            (*"EqAllAllSimple" -> EqAllAllSimple,*)
            "RulesEntryIn"-> RulesEntryIn,
            "RulesExitValues" -> RulesExitValues,
            (*"EqAllAllRules" -> EqAllAllRules,*)
            "Nlhs" -> Nlhs,
            "nocasesystem" -> nocasesystem, 
            "nocaserules" -> nocaserules,
            "MinimalTimeRhs" -> MinimalTimeRhs,
            "EqCriticalCase" -> EqCriticalCase ,
            "criticalreduced1" -> criticalreduced1,
            (*"criticalreduced2" -> criticalreduced2,*)
            "Nrhs" -> Nrhs,
            "EqGeneralCase" -> EqGeneralCase,
            "TOL"->TOL
            }],
            Association[{
                "Data" -> Data,
              (*Graph structure*)
              (*"BG" -> BG, 
              "EntranceVertices" -> EntranceVertices, 
              "InwardVertices" -> InwardVertices, 
              "InEdges" -> InEdges, 
              "ExitVertices" -> ExitVertices, 
              "OutwardVertices" -> OutwardVertices, 
              "OutEdges" -> OutEdges, 
              "AuxiliaryGraph" -> AuxiliaryGraph,*) 
              "FG" -> FG, (*
              "VL" -> VL, 
              "EL" -> EL, 
              "BEL" -> BEL, 
              "FVL" -> FVL, 
              "AllTransitions" -> AllTransitions, 
              "NoDeadEnds" -> NoDeadEnds, 
              "NoDeadStarts" -> NoDeadStarts, 
              (*variables*)
              "jargs" -> jargs, 
              "js" -> js, 
              "jvars" -> jvars, 
              "jts" -> jts, 
              "jtvars" -> jtvars, 
              "uargs" -> uargs, 
              "us" -> us, 
              "uvars" -> uvars,
              "SwitchingCosts" -> SwitchingCosts, 
              "OutRules" -> OutRules, 
              "InRules" -> InRules, 
              
              
              
              "EntryArgs" -> EntryArgs, 
              "EntryDataAssociation" -> EntryDataAssociation, 
              "jays" -> SignedCurrents, 
              "NonZeroEntryCurrents" -> NonZeroEntryCurrents, 
              "ExitCosts" -> ExitCosts, 
              
              
              *) 
              (*equations*)
              (*complementarity*)
            (*          "EqCurrentCompCon" -> EqCurrentCompCon, 
              "EqTransitionCompCon" -> EqTransitionCompCon, 
              "EqCompCon" -> EqCompCon, 
              "EqAllComp" -> EqAllComp, (*union of all complementarity conditions*)*)
              
              (*linear equations (and inequalities)*)
              (*"EqPosCon" -> EqPosCon, 
              "EqBalanceSplittingCurrents" -> EqBalanceSplittingCurrents, 
              "EqBalanceGatheringCurrents" -> EqBalanceGatheringCurrents, 
              "EqEntryIn" -> EqEntryIn, 
              "EqExitValues" -> EqExitValues, 
              "EqSwitchingConditions" -> EqSwitchingConditions, 
              "EqValueAuxiliaryEdges" -> EqValueAuxiliaryEdges, 
              "EqAll" -> EqAll, *)
             (* "EqAllCompRules" -> EqAllCompRules,
              "EqAllRules" -> EqAllRules,*)
              "EqAllAll" -> EqAllAll,
              "BoundaryRules" -> BoundaryRules,
              "reduced1" -> reduced1,
              "reduced2" -> reduced2,
            (*          "EqAllAllSimple" -> EqAllAllSimple,*)
            (*          "RulesEntryIn"-> RulesEntryIn,
              "RulesExitValues" -> RulesExitValues,*)
            (*          "EqAllAllRules" -> EqAllAllRules,*)
              (*"reduced system" -> reduced[[1]],
              "reducing rules" -> reduced[[2]],*)
              "Nlhs" -> Nlhs,
              "EqCriticalCase" -> EqCriticalCase ,
              "criticalreduced1" -> criticalreduced1,
              "criticalreduced2" -> criticalreduced2,
              "Nrhs" -> Nrhs,
              "EqGeneralCase" -> EqGeneralCase,
              "TOL" -> TOL
              }]
        ]
    ]

CriticalBundle[Data_Association] :=
    Module[ {d2e},
        d2e = D2E[Data];
        CriticalCongestionSolver[d2e]
    ]


CriticalCongestionSolverN[D2E_Association] :=
    Module[ {system, rules, sol,
        EqAllAll = Lookup[D2E, "EqAllAll", Print["No equations to solve."];
                                           Return[]], 
        EqCriticalCase = Lookup[D2E,"EqCriticalCase", Print["Critical case equations are missing."];
                                                      Return[]]
        },
        rules = First @ Solve[EqCriticalCase];
        system = EqAllAll /. rules;
        {system, rules} = CleanEqualities[{system, rules}];
        {system, rules} = {system, rules}/.rules;
        system = NewReduce[BooleanConvert[system,"CNF"]];
        {system, rules} = CleanEqualities[{system, rules}];
        rules = Simplify @ rules;
        Print[system];
        Which[system === True || Head[system] === Equal || (Head[system] === And && DeleteDuplicates[Head /@ system] === Equal),
            sol = Solve[system];
            If[ Length[sol] == 1,
                {system, rules} = {system, AssociateTo[ rules, First@ sol]} /. First@ sol,
                Print[sol];
            ],
            system === False,
            Print["CriticalCongestionSolver: ",Style["Incompatible system", Red],"."],
            True,
            Print["CriticalCongestionSolver: The system does not have the original structure: \n", system];
            {system,rules} = CleanEqualities[{system, rules}];
            system = Reduce[system];
            Print["CriticalCongestionSolver: ", Style["Multiple solutions", Red]," \n\t\t", system(*, "\n\tFinding one instance..."*)];
            (*FindInstance!!!*)
            (*Module[ {solved = And @@ (Outer[Equal,Keys[rules],Values[rules]]//Diagonal) , onesol},
                onesol = First@FindInstance[system && solved, Join[js, us, jts], Reals];    
                (*TODO: should we choose a method to select a solution?*)
                rules = Association[onesol];
                system = system/.rules
            ];*)
        ];
        If[ system === True,
            Print["CriticalCongestionSolver: Critical congestion solved."];,
            Print["CriticalCongestionSolver: There are ", Style["multiple solutions", Red]," for the given data."];
        ];
        {system, Simplify /@ rules}
    ]
    
CriticalCongestionSolver[D2E_Association] :=
    Module[ {system, rules,
        EqAllAll = Lookup[D2E, "EqAllAll", Print["No equations to solve."];
                                           Return[]], 
        EqCriticalCase = Lookup[D2E,"EqCriticalCase", Print["Critical case equations are missing."];
                                                      Return[]],
        InitRules = Lookup[D2E, "BoundaryRules", {}]                                              
        },
        EqAllAll = EqAllAll /. InitRules;
        EqCriticalCase = EqCriticalCase /. InitRules;
        rules = First @ Solve @ EqCriticalCase//Quiet;
        system = EqAllAll/.rules;
        {system,rules} = FixedPoint[CriticalCongestionStep, {system,Join[InitRules,rules]}];
        rules = Simplify /@ rules;
        Simplify /@ ({system, rules}/.rules)
        
    ];
    
    CriticalCongestionStep[{sys_Or, rul_}]:= 
    	{Reduce[sys, Reals], rul}
    
    
    CriticalCongestionStep[{False, rul_}] :=
        {False, rul}
    
    CriticalCongestionStep[{sys_,rul_}] :=
        Module[ {system = sys,rules = rul},
            {system,rules} = CleanEqualities[{system, rules}];
            system = BooleanConvert[system,"CNF"];
            system = NewReduce[system];
            {system,rules}
        ]
    
End[]