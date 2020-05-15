(* Wolfram Language package *)
(*Assembles the stationary mfg equations using the notion of current*)
DataToEquations::usage = "takes data assotiation to the equations";

Begin["`Private`"]
(*TODO remove (most) functions from the association DONE*)
(*TODO how do we use usage message to make template on the notebook*)

DataToEquations::usage = 
"DataToEquations[<|\"Vertices List\"->{x}|>] Inputs the data as an association and returns the equations for the stationary mean-field game on the network, except for the coupling equations.";

DataToEquations[Data_?AssociationQ] = 
 Module[{BG, EntranceVertices, InwardVertices, InEdges, ExitVertices, 
   OutwardVertices, OutEdges, AuxiliaryGraph, FG, VL, EL, BEL, jargs, 
   js, jvars, FVL, AllTransitions, jts, jtvars, uargs, us, uvars, 
   EqPosCon, CurrentCompCon, EqCurrentCompCon, TransitionCompCon, 
   EqTransitionCompCon, IncomingEdges, OutgoingEdges, NoDeadEnds, 
   NoDeadStarts, CurrentSplitting, CurrentGathering, 
   EqBalanceSplittingCurrents, EqBalanceGatheringCurrents, EntryArgs, 
   EntryDataAssociation, EqEntryIn, NonZeroEntryCurrents, ExitCosts, 
   ExitValues, EqExitValues, SwitchingCosts, OutRules, InRules, 
   Transu, EqSwitchingConditions, Compu, EqCompCon, EqValue, 
   EqValueAuxiliaryEdges, EqAllComp, EqAll},
  BG = AdjacencyGraph[Data["Vertices List"], Data["Adjacency Matrix"],
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
  AllTransitions = 
   TransitionsAt[FG, #] & /@ FVL // 
    Catenate(*at vertex from first edge to second edge*);
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
  (*EqAll=EqPosCon&&EqCurrentCompCon&&EqTransitionCompCon;*)
  
  IncomingEdges[k_] := {k, #1} & /@ 
    IncidenceList[FG, 
     k](*all edges "oriented" towards the vertex k*)
  (*all edges\ \
"oriented" away from the vertex k*);
  OutgoingEdges[k_] := 
   OtherWay /@ ({k, #} & /@ IncidenceList[FG, k]);
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
  (*EqAll=EqAll&&EqBalanceSplittingCurrents&&
  EqBalanceGatheringCurrents;*)
  
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
  Compu[{v_, edge1_, edge2_}] := (jtvars[{v, edge1, edge2}] == 0) || 
    uvars[{v, edge1}] - uvars[{v, edge2}] - 
      SwitchingCosts[{v, edge1, edge2}] == 0;
  EqCompCon = And @@ Compu /@ AllTransitions;
  (*EqAll=EqAll&&EqCompCon;*)
  (*EqValue=
  And@@((DifferenceU[#]\[Equal]uvars[AtHead[#]]-uvars[AtTail[#]])&/@
  EdgeList[BG]);*)(*use the same reasoning to approach the more \
general cases.*)
  EqValueAuxiliaryEdges = 
   And @@ ((uvars[AtHead[#]] - uvars[AtTail[#]] == 0) & /@ 
      EdgeList[AuxiliaryGraph]);
  (*EqNoInt=EqAll&&EqValueAuxiliaryEdges;*)
  
  EqAllComp = EqCurrentCompCon && EqTransitionCompCon && EqCompCon;
  EqAll = EqAll && EqValueAuxiliaryEdges;
  (*EqAll=EqAll&&EqValue&&EqValueAuxiliaryEdges;*)
  
  Association[{
  	"BG" -> BG, "EntranceVertices" -> EntranceVertices, 
    "InwardVertices" -> InwardVertices, "InEdges" -> InEdges, 
    "ExitVertices" -> ExitVertices, 
    "OutwardVertices" -> OutwardVertices, "OutEdges" -> OutEdges, 
    "AuxiliaryGraph" -> AuxiliaryGraph, "FG" -> FG, "VL" -> VL, 
    "EL" -> EL, "BEL" -> BEL, "jargs" -> jargs, "js" -> js, 
    "jvars" -> jvars, "FVL" -> FVL, 
    "AllTransitions" -> AllTransitions, "jts" -> jts, 
    "jtvars" -> jtvars, "uargs" -> uargs, "us" -> us, 
    "uvars" -> uvars, "EqPosCon" -> EqPosCon, 
    "CurrentCompCon[a\[DirectedEdge]b]" -> CurrentCompCon, 
    "EqCurrentCompCon" -> EqCurrentCompCon, 
    "TransitionCompCon[{v, edge1, edge2}]" -> TransitionCompCon, 
    "EqTransitionCompCon" -> EqTransitionCompCon, 
    "IncomingEdges[k]" -> IncomingEdges, 
    "OutgoingEdges[k]" -> OutgoingEdges, "NoDeadEnds" -> NoDeadEnds, 
    "NoDeadStarts" -> NoDeadStarts, 
    "CurrentSplitting[{c,a\[DirectedEdge]b}]" -> CurrentSplitting, 
    "CurrentGathering[{c,a\[DirectedEdge]b}]" -> CurrentGathering, 
    "EqBalanceSplittingCurrents" -> EqBalanceSplittingCurrents, 
    "EqBalanceGatheringCurrents" -> EqBalanceGatheringCurrents, 
    "EntryArgs" -> EntryArgs, 
    "EntryDataAssociation" -> EntryDataAssociation, 
    "EqEntryIn" -> EqEntryIn, 
    "NonZeroEntryCurrents" -> NonZeroEntryCurrents, 
    "ExitCosts" -> ExitCosts, 
    "ExitValues[a\[DirectedEdge]b]" -> ExitValues, 
    "EqExitValues" -> EqExitValues, 
    "SwitchingCosts" -> SwitchingCosts, "OutRules" -> OutRules, 
    "InRules" -> InRules, "Transu[{v,edge1,edge2}]" -> Transu, 
    "EqSwitchingConditions" -> EqSwitchingConditions, 
    "Compu[{v,edge1,edge2}]" -> Compu, "EqCompCon" -> EqCompCon, 
    "EqValue" -> EqValue, 
    "EqValueAuxiliaryEdges" -> EqValueAuxiliaryEdges, 
    "EqAllComp" -> EqAllComp, "EqAll" -> EqAll
    
    }]
    ]

(*DataToEquations[Data_?AssociationQ] :=
Module[{
	GraphData = DataToGraph[Data], 
    VL  = Data["Vertices List"], 
    EL, 
    BEL, 
    jargs, 
    js, 
    jvars, 
    FVL, 
    AllTransitions, 
    jts, 
    jtvars, 
    uargs, 
    us, 
    uvars, 
    FG, 
    BG, 
    EqPosCon},
    BG = GraphData["BG"];
     EL = EdgeList[FG];
     BEL = EdgeList[BG];
     jargs = Join[AtHead /@ EL, AtTail /@ EL];
     js = Unique["j"] & /@ jargs;
     jvars = AssociationThread[jargs, js];
     FVL = VertexList[FG];
     AllTransitions = 
      TransitionsAt[FG, #] & /@ FVL // 
    Catenate(*at vertex from first edge to second edge*);
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
     EqAll = EqPosCon && EqCurrentCompCon && EqTransitionCompCon;
     IncomingEdges[k_] :=
         {k, #1} & /@ 
         IncidenceList[FG,k](*all edges "oriented" towards the vertex k*)
     (*all edges "oriented" away from the vertex k*);
     OutgoingEdges[k_] :=
         OtherWay /@ ({k, #} & /@ IncidenceList[FG, k]);
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
         NoDeadStarts );
     EqAll = 
      EqAll && EqBalanceSplittingCurrents && EqBalanceGatheringCurrents;
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
         Total[ uvars /@ {{b, DirectedEdge[a, b]}}] == ExitCosts[b];
     EqExitValues = 
      And @@ (ExitValues /@ 
         IncidenceList[AuxiliaryGraph, 
          OutwardVertices /@ ExitVertices]);
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
     Transu[{v_, edge1_, edge2_}] :=
         NonNegative[
          SwitchingCosts[{v, edge2, edge1}] + uvars[{v, edge1}] - 
           uvars[{v, edge2}]];
     EqSwitchingConditions = (And @@ Transu /@ AllTransitions);
     EqAll = EqAll && EqSwitchingConditions;
     Compu[{v_, edge1_, edge2_}] :=
         (jtvars[{v, edge1, edge2}] == 0) || 
         uvars[{v, edge1}] - uvars[{v, edge2}] - 
         SwitchingCosts[{v, edge1, edge2}] == 0;
     EqCompCon = And @@ Compu /@ AllTransitions;
     EqAll = EqAll && EqCompCon;
   
     EqValue = 
      And @@ ((DifferenceU[#] == uvars[AtHead[#]] - uvars[AtTail[#]]) & /@
          EdgeList[BG]);(*use the same reasoning to approach the more general \
   cases.*)
     EqValueAuxiliaryEdges = 
      And @@ ((uvars[AtHead[#]] - uvars[AtTail[#]] == 0) & /@ 
         EdgeList[AuxiliaryGraph]);
     EqNoInt = EqAll && EqValueAuxiliaryEdges;
     EqAll = EqAll && EqValue && EqValueAuxiliaryEdges
    ];*)


(*TODO is it possible to assamble the equations by edge/vertex? Then we could try to difuse the information through the graph.*)
(*EqAll=DataToEquations[Data];*)

End[]

