(* Wolfram Language package *)

(*Assembles the stationary mfg equations using the notion of current*)
DataToEquations::usage = 
"DataToEquations[<|\"Vertices List\" -> {1, 2, 3}, \"Adjacency Matrix\" -> {{0, 0, 0}, {1, 0, 1}, {0, 0, 0}}, 
 \"Entrance Vertices and Currents\" -> {{2, I1}}, 
 \"Exit Vertices and Terminal Costs\" -> {{1, U1}, {3, U2}}, 
 \"Switching Costs\" -> {{1, 2, 3, S1}, {3, 2, 1, S2}}|>] returns the equations for the stationary mean-field game on the network, except for the coupling equations.";

Begin["`Private`"]



DataToEquations[Data_?AssociationQ] :=
    Module[ {BG, EntranceVertices, InwardVertices, InEdges, ExitVertices, 
      OutwardVertices, OutEdges, AuxiliaryGraph, FG, VL, EL, BEL, jargs, 
      js, jvars, FVL, AllTransitions, jts, jtvars, uargs, us, uvars, 
      EqPosCon, CurrentCompCon, EqCurrentCompCon, TransitionCompCon, 
      EqTransitionCompCon, IncomingEdges, OutgoingEdges, NoDeadEnds, 
      NoDeadStarts, CurrentSplitting, CurrentGathering, 
      EqBalanceSplittingCurrents, EqBalanceGatheringCurrents, EntryArgs, 
      EntryDataAssociation, EqEntryIn, NonZeroEntryCurrents, ExitCosts, 
      ExitValues, EqExitValues, SwitchingCosts, OutRules, InRules, 
      Transu, EqSwitchingConditions, Compu, EqCompCon,
      EqValueAuxiliaryEdges, EqAllComp, EqAll, SignedCurrents, Nrhs, Nlhs, RulesEntryIn, 
      RulesExitValues, EqAllRules, EqAllCompRules, EqAllAll, EqAllAllSimple, EqAllAllRules, reduced,
      EqCriticalCase, EqCriticalCaseRules, criticalreduced },
        BG = AdjacencyGraph[
            Data["Vertices List"], 
            Data["Adjacency Matrix"], 
            VertexLabels -> "Name"];
        EntranceVertices = First /@ Data["Entrance Vertices and Currents"];
        InwardVertices = 
         AssociationThread[EntranceVertices, 
          Unique["en"] & /@ EntranceVertices];
        (*Defines auxiliary vertices for the entrance vertices*)
        InEdges = 
         MapThread[
          DirectedEdge, {InwardVertices /@ EntranceVertices, 
           EntranceVertices}];
        ExitVertices = First /@ Data["Exit Vertices and Terminal Costs"];
        OutwardVertices = 
         AssociationThread[ExitVertices, Unique["ex"] & /@ ExitVertices];
        OutEdges = 
         MapThread[
          DirectedEdge, {ExitVertices, OutwardVertices /@ ExitVertices}];
        AuxiliaryGraph = 
         Graph[Join[InEdges, OutEdges], VertexLabels -> "Name"];
        FG = EdgeAdd[BG, Join[InEdges, OutEdges]];
        VL = Data["Vertices List"];
        EL = EdgeList[FG];
        BEL = EdgeList[BG];
        jargs = Join[AtHead /@ EL, AtTail /@ EL];
        js = Unique["j"] & /@ jargs;
        jvars = AssociationThread[jargs, js];
        FVL = VertexList[FG];
        SignedCurrents =   AssociationThread[BEL, (jvars[AtHead[#]] - jvars[AtTail[#]] &) /@ BEL];
        AllTransitions = 
         TransitionsAt[FG, #] & /@ FVL // Catenate(*at vertex from first edge to second edge*);
        jts = Unique["jt"] & /@ AllTransitions;
        jtvars = AssociationThread[AllTransitions, jts];
        uargs = Join[AtTail /@ EL, AtHead /@ EL];
        us = Unique["u"] & /@ uargs;
        uvars = AssociationThread[uargs, us];
        EqPosCon = And @@ (NonNegative /@ Join[jvars, jtvars]);
        CurrentCompCon[a_ \[DirectedEdge] b_] :=
            jvars[{a, a \[DirectedEdge] b}] == 0 || 
             jvars[{b, a \[DirectedEdge] b}] == 0;
        EqCurrentCompCon = And @@ (CurrentCompCon /@ EL);
        TransitionCompCon[{v_, edge1_, edge2_}] :=
            jtvars[{v, edge1, edge2}] == 0 || jtvars[{v, edge2, edge1}] == 0; 
        EqTransitionCompCon = 
         And @@ ((Sort /@ TransitionCompCon /@ AllTransitions) // Union);
 
        NoDeadEnds = IncomingEdges /@ VL // Flatten[#, 1] &;
        NoDeadStarts = OutgoingEdges /@ VL // Flatten[#, 1] &;
        CurrentSplitting[{c_, a_ \[DirectedEdge] b_}] :=
            Select[AllTransitions, (Take[#, 2] == {c, 
                 a \[DirectedEdge] b}) &];
        CurrentGathering[{c_, a_ \[DirectedEdge] b_}] :=
            Select[AllTransitions, (Part[#, {1, 3}] == 
                OtherWay[{c, a \[DirectedEdge] b}]) &];
        EqBalanceSplittingCurrents = 
         And @@ ((jvars[#] == Total[jtvars /@ CurrentSplitting[#]]) & /@ 
            NoDeadEnds);
        EqBalanceGatheringCurrents = 
         And @@ ((jvars[#] == Total[jtvars /@ CurrentGathering[#]]) & /@ 
            NoDeadStarts);
        EqAll = 
         EqPosCon && EqBalanceSplittingCurrents && 
          EqBalanceGatheringCurrents;
   
        EntryArgs = 
         AtHead /@ ((EdgeList[
                AuxiliaryGraph, _ \[DirectedEdge] #] & /@ (First /@ 
                Data["Entrance Vertices and Currents"])) // Flatten[#, 1] &);
        EntryDataAssociation = 
         AssociationThread[EntryArgs, 
          Last /@ Data["Entrance Vertices and Currents"]];
        EqEntryIn = 
         And @@ ((jvars[#] == EntryDataAssociation[#]) & /@ (AtHead /@ 
              InEdges));
        RulesEntryIn = ToRules[EqEntryIn];
        NonZeroEntryCurrents = 
         And @@ (Positive[EntryDataAssociation[#]] & /@ (AtHead /@ InEdges));
        (*useful for the general case...*)
        EqAll = EqAll && EqEntryIn;
        ExitCosts = 
         AssociationThread[
          OutwardVertices /@ (First /@ 
             Data["Exit Vertices and Terminal Costs"]), 
          Last /@ Data["Exit Vertices and Terminal Costs"]];
        ExitValues[a_ \[DirectedEdge] b_] :=
            Total[uvars /@ {{b, DirectedEdge[a, b]}}] == ExitCosts[b];
        EqExitValues = 
         And @@ (ExitValues /@ 
            IncidenceList[AuxiliaryGraph, 
             OutwardVertices /@ ExitVertices]);
        RulesExitValues = ToRules[EqExitValues];
        EqAll = EqAll && EqExitValues;
        SwitchingCosts = 
         AssociationThread[
          Join[AllTransitions, 
           triple2path[Take[#, 3], FG] & /@ Data["Switching Costs"]], 
          Join[Table[0, Length[AllTransitions]], 
           Last[#] & /@ 
            Data["Switching Costs"]]](*Swithing cost is initialized with 0. \
      AssociationThread associates the last association!*);
        OutRules = 
         Rule[#, Infinity] & /@ (Outer[Flatten[{AtTail[#1], #2}] &, 
             OutEdges, EL] // Flatten[#, 1] &);
        AssociateTo[SwitchingCosts, Association[OutRules]];
        InRules = 
         Rule[#, Infinity] & /@ (Outer[{#2[[2]], #1, #2} &, 
             IncidenceList[FG, #] & /@ EntranceVertices, InEdges] // 
      Flatten[#, 2] &);
        AssociateTo[SwitchingCosts, Association[InRules]];
        Transu[{v_, edge1_, edge2_}] :=
            uvars[{v, edge2}] <= 
             SwitchingCosts[{v, edge2, edge1}] + uvars[{v, edge1}];
        (*Transu[{v_,edge1_,edge2_}]:=NonNegative[SwitchingCosts[{v,edge2,
        edge1}]+uvars[{v,edge1}]-uvars[{v,edge2}]];*)
        EqSwitchingConditions = (And @@ Transu /@ AllTransitions);
        EqAll = EqAll && EqSwitchingConditions;
        Compu[{v_, edge1_, edge2_}] :=
            (jtvars[{v, edge1, edge2}] == 0) || 
            uvars[{v, edge1}] - uvars[{v, edge2}] - 
            SwitchingCosts[{v, edge1, edge2}] == 0;
        EqCompCon = And @@ Compu /@ AllTransitions;
        EqValueAuxiliaryEdges = And @@ ((uvars[AtHead[#]] - uvars[AtTail[#]] == 0) & /@ EdgeList[AuxiliaryGraph]);

        EqAllComp = EqCurrentCompCon && EqTransitionCompCon && EqCompCon;
        EqAll = EqAll && EqValueAuxiliaryEdges;
        
        (*Pre-processing: substitute rules in all equations and hand-sides...*)
        EqAllCompRules = EqAllComp /. Join[RulesEntryIn, RulesExitValues];
        EqAllRules = EqAll /. Join[RulesEntryIn, RulesExitValues];
        EqAllAll = EqAll && EqAllComp;(*this is in the place of Boo, without the brackets.*)
        EqAllAllRules = EqAllCompRules && EqAllRules;
       
        reduced = 
  			FixedPoint[EqEliminatorX, {EqAllAll, Join[RulesEntryIn, RulesExitValues]}, 10];
 		EqAllAllSimple = reduced[[1]];
 		
        Nlhs = Flatten[uvars[AtHead[#]] - uvars[AtTail[#]] + SignedCurrents[#] & /@ BEL];
        Nrhs =  Flatten[Intg[SignedCurrents[#]] + SignedCurrents[#] & /@ BEL];
       	Print["more stuff"];
       	EqCriticalCase = # && And @@ (# == 0) & /@ Nlhs;
       	EqCriticalCaseRules = EqCriticalCase /. reduced[[2]];
       
       	
       	criticalreduced = FixedPoint[EqEliminatorX, {reduced[[1]] && EqCriticalCaseRules, reduced[[2]]},10]; 
       	(*TODO see if this is enough in the place of startsolverX*)
       	
        Association[{
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
          (*equations*)
          (*complementarity*)
          "CurrentCompCon" -> CurrentCompCon, 
          "EqCurrentCompCon" -> EqCurrentCompCon, 
          "TransitionCompCon" -> TransitionCompCon, 
          "EqTransitionCompCon" -> EqTransitionCompCon, 
          "Compu" -> Compu, 
          "EqCompCon" -> EqCompCon, 
          "EqAllComp" -> EqAllComp, (*union of all complementarity conditions*)
          (*linear equations (and inequalities)*)
          "EqPosCon" -> EqPosCon, 
          "CurrentSplitting" -> CurrentSplitting, 
          "CurrentGathering" -> CurrentGathering, 
          "EqBalanceSplittingCurrents" -> EqBalanceSplittingCurrents, 
          "EqBalanceGatheringCurrents" -> EqBalanceGatheringCurrents, 
          "EntryArgs" -> EntryArgs, 
          "EntryDataAssociation" -> EntryDataAssociation, 
          "EqEntryIn" -> EqEntryIn, 
          "NonZeroEntryCurrents" -> NonZeroEntryCurrents, 
          "ExitCosts" -> ExitCosts, 
          "ExitValues" -> ExitValues, 
          "EqExitValues" -> EqExitValues, 
          "SwitchingCosts" -> SwitchingCosts, 
          "OutRules" -> OutRules, 
          "InRules" -> InRules, 
          "Transu" -> Transu, 
          "EqSwitchingConditions" -> EqSwitchingConditions, 
          "EqValueAuxiliaryEdges" -> EqValueAuxiliaryEdges, 
          "EqAll" -> EqAll, 
          "jays" -> SignedCurrents, 
          "Nrhs" -> Nrhs, 
          "Nlhs" -> Nlhs,
          "EqAllCompRules" -> EqAllCompRules,
          "EqAllRules" -> EqAllRules,
          "EqAllAll" -> EqAllAll,
          "EqAllAllSimple" -> EqAllAllSimple,
          "RulesEntryIn"-> RulesEntryIn,
          "RulesExitValues" -> RulesExitValues,
          "EqAllAllRules" -> EqAllAllRules,
          "reduced system" -> reduced[[1]],
          "reducing rules" -> reduced[[2]],
          "EqCriticalCaseRules" -> EqCriticalCaseRules,
          "EqCriticalCase" -> EqCriticalCase 
          }]
    ]

End[]

