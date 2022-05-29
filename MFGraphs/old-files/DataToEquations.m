(*Wolfram Language package*)
D2E::usage = "D2E[<|\"Vertices List\" -> {1, 2, 3}, \"Adjacency \
Matrix\" -> {{0, 0, 0}, {1, 0, 1}, {0, 0, 0}}, 
 \"Entrance Vertices and Currents\" -> {{2, I1}}, 
 \"Exit Vertices and Terminal Costs\" -> {{1, U1}, {3, U2}}, 
 \"Switching Costs\" -> {{1, 2, 3, S1}, {3, 2, 1, S2}}|> returns the \
equations for the stationary mean-field game on the network.]"

Begin["`Private`"]
D2E[Data_Association] := 
 Module[{BG, FG, AuxiliaryGraph, VL, FVL, EntranceVertices, 
   InwardVertices, ExitVertices, OutwardVertices, EL, BEL, OutEdges, 
   InEdges, ExitNeighbors, AllTransitions, EqCcs, 
   SC = Lookup[Data, "Switching Costs", {}], jargs, js, jvars, jts, 
   jtvars, uargs, us, uvars, SignedCurrents, EntryArgs, NoDeadEnds, 
   NoDeadStarts, ExitCosts, EntryDataAssociation, SwitchingCosts, 
   EqPosCon, EqCurrentCompCon, EqTransitionCompCon, 
   EqSwitchingByVertex, EqSwitchingConditions, EqCompCon, 
   EqValueAuxiliaryEdges, EqBalanceSplittingCurrents, 
   EqBalanceGatheringCurrents, Nrhs, Nlhs, MinimalTimeRhs, AllOr, 
   EqAllAll, AllIneq, EqCriticalCase, EqMinimalTime, EqNonCritical, 
   EqPosJs, TrueEq, costpluscurrents, RuleBalanceGatheringCurrents, 
   RuleEntryIn, RuleEntryOut, RuleExitValues, RuleExitCurrentsIn, 
   InitRules, OutRules, InRules, RuleNonCritical, RuleNonCritical1, 
   RulesCriticalCase, RulesCriticalCase1,(*,EqNonCritical1*)EqPosJts, 
   BalanceSplittingCurrents, BalanceGatheringCurrents, EqEntryIn, 
   Kirchhoff, RuleValueAuxiliaryEdges, BM, KM, 
   vars},(*Checking consistency on the swithing costs*)
  If[SC =!= {}, 
   ConsistentSwithingCosts[sc_][{a_, b_, c_, S_}] := 
    Module[{or, de, bounds}, or = Cases[sc, {a, b, _, _}];
     de = Cases[sc, {_, b, c, _}];
     bounds = Outer[Plus, Last /@ or, Last /@ de] // Flatten;
     (*And@@(NonNegative[#-S]&)/@bounds*)
     And @@ (# >= S &) /@ bounds];
   EqCcs = ConsistentSwithingCosts[SC] /@ SC;
   EqCcs = Simplify[EqCcs];
   If[AnyTrue[EqCcs, # === False &], 
    Print[
     "DataToEquations: Triangle inequalities for switching costs: \n",
      Select[AssociationThread[SC, EqCcs], 
      Function[exp, ! TrueQ[exp]]]];
    Print["DataToEquations: The switching costs are ", 
     Style["incompatible", Red], ". \nStopping!"];
    Return[]]];
  (****Graph stuff****)
  BG = AdjacencyGraph[Data["Vertices List"], Data["Adjacency Matrix"],
     VertexLabels -> "Name", DirectedEdges -> True];
  EntranceVertices = First /@ Data["Entrance Vertices and Currents"];
  ExitVertices = First /@ Data["Exit Vertices and Terminal Costs"];
  Clear["en*"];
  (*InwardVertices defines auxiliary vertices for the entrance \
vertices*)
  InwardVertices = 
   AssociationThread[EntranceVertices, 
    Symbol["en" <> ToString[#]] & /@ EntranceVertices];
  Clear["ex*"];
  OutwardVertices = 
   AssociationThread[ExitVertices, 
    Symbol["ex" <> ToString[#]] & /@ ExitVertices];
  (*InEdges defines auxiliary arguments for the entrance vertices*)
  InEdges = 
   MapThread[
    DirectedEdge, {InwardVertices /@ EntranceVertices, 
     EntranceVertices}];
  OutEdges = 
   MapThread[
    DirectedEdge, {ExitVertices, OutwardVertices /@ ExitVertices}];
  AuxiliaryGraph = 
   Graph[Join[InEdges, OutEdges], VertexLabels -> "Name", 
    GraphLayout -> "SpringEmbedding"];
  FG = EdgeAdd[BG, Join[InEdges, OutEdges]];
  VL = Data["Vertices List"];
  EL = EdgeList[FG];
  BEL = EdgeList[BG];
  FVL = VertexList[FG];
  (*arguments*)
  jargs = Flatten[#, 1] &@({AtTail@#, AtHead@#} & /@ EL);
  uargs = jargs;
  AllTransitions = 
   TransitionsAt[FG, #] & /@ FVL // 
    Catenate(*at vertex from first edge to second edge*);
  EntryArgs = 
   AtHead /@ ((EdgeList[
          AuxiliaryGraph, _ \[DirectedEdge] #] & /@ (First /@ 
          Data["Entrance Vertices and Currents"])) // Flatten[#, 1] &);
  EntryDataAssociation = 
   RoundValues@
    AssociationThread[EntryArgs, 
     Last /@ Data["Entrance Vertices and Currents"]];
  ExitCosts = 
   AssociationThread[
    OutwardVertices /@ (First /@ 
       Data["Exit Vertices and Terminal Costs"]), 
    Last /@ Data["Exit Vertices and Terminal Costs"]];
  (*variables*)
  js = Table[Symbol["j" <> ToString[k]], {k, 1, Length@jargs}];
  jvars = AssociationThread[jargs, js];
  jts = Table[
    Symbol["jt" <> ToString[k]], {k, 1, Length@AllTransitions}];
  jtvars = AssociationThread[AllTransitions, jts];
  us = Table[Symbol["u" <> ToString[k]], {k, 1, Length@uargs}];
  uvars = AssociationThread[uargs, us];
  costpluscurrents = 
   Table[Symbol["cpc" <> ToString[k]], {k, 1, Length@BEL}];
  (*costpluscurrentsvars=AssociationThread[BEL,costpluscurrents];*)
  SignedCurrents = 
   AssociationThread[
    BEL, (jvars[AtHead[#]] - jvars[AtTail[#]] &) /@ BEL];
  Print["D2E: Variables are all set"];
  (*Elements of the system*)(*Swithing cost is initialized with 0. \
AssociationThread associates the last association!*)
  SwitchingCosts = 
   AssociationThread[
    Join[AllTransitions, triple2path[Take[#, 3], FG] & /@ SC], 
    Join[0 & /@ AllTransitions, Last[#] & /@ SC]];
  SC = path2triple /@ Normal[SwitchingCosts];
  If[SC =!= {}, Print["Second check:"];
   (*ConsistentSwithingCosts[{a_,b_,c_,S_}]:=Module[{sc=SC,or,de,
   bounds},or=Cases[sc,{a,b,_,_}];
   de=Cases[sc,{_,b,c,_}];
   bounds=Outer[Plus,Last/@or,Last/@de]//Flatten;
   (*And@@(NonNegative[#-S]&)/@bounds*)And@@(#>=S&)/@bounds];*)
   EqCcs = ConsistentSwithingCosts[SC] /@ SC;
   EqCcs = Simplify[EqCcs];
   (*Print[EqCcs];*)
   If[AnyTrue[EqCcs, # === False &], 
    Print[
     "DataToEquations: Triangle inequalities for switching costs: ", 
     Select[AssociationThread[SC, EqCcs], 
      Function[exp, ! TrueQ[exp]]]];
    Print["DataToEquations: The switching costs are ", 
     Style["incompatible", Red], ". \nStopping!"];
    Return[]]];
  EqPosJs = And @@ (# >= 0 & /@ Join[jvars]);(*Inequality*)
  EqPosJts = And @@ (# >= 0 & /@ Join[jtvars]);(*Inequality*)
  EqPosCon = EqPosJs && EqPosJts;(*Inequality*)
  EqCurrentCompCon = And @@ (CurrentCompCon[jvars] /@ EL);(*Or*)
  EqTransitionCompCon = 
   And @@ ((Sort /@ TransitionCompCon[jtvars] /@ AllTransitions) // 
      Union);(*Or*)(*Balance Splitting Currents in the full graph*)
  NoDeadEnds = IncomingEdges[FG] /@ VL // Flatten[#, 1] &;
  EqBalanceSplittingCurrents = 
   And @@ ((jvars[#] == 
         Total[jtvars /@ CurrentSplitting[AllTransitions][#]]) & /@ 
      NoDeadEnds);(*Equal*)
  BalanceSplittingCurrents = ((jvars[#] - 
        Total[jtvars /@ CurrentSplitting[AllTransitions][#]]) & /@ 
     NoDeadEnds);
  (*Gathering currents in the inside of the basic graph*)
  NoDeadStarts = OutgoingEdges[FG] /@ VL // Flatten[#, 1] &;
  RuleBalanceGatheringCurrents = (jvars[#] -> 
       Total[jtvars /@ CurrentGathering[AllTransitions][#]]) & /@ 
    NoDeadStarts;(*Rule*)(*First rules:
  these have some j in terms of jts*)
  InitRules = Association[RuleBalanceGatheringCurrents];
  (*get equations for the exit currents at the entry vertices*)
  BalanceGatheringCurrents = ((-jvars[#] + 
        Total[jtvars /@ CurrentGathering[AllTransitions][#]]) & /@ 
     NoDeadStarts);
  EqBalanceGatheringCurrents = 
   Simplify /@ (And @@ (# == 0 & /@ BalanceGatheringCurrents));
  (*Incoming currents*)
  EqEntryIn = (jvars[#] == EntryDataAssociation[#]) & /@ (AtHead /@ 
      InEdges);(*List of Equals*)
  RuleEntryIn = Flatten[ToRules /@ EqEntryIn];(*List of Rules*)
  Kirchhoff = 
   Join[EqEntryIn, (# == 0 & /@ (BalanceGatheringCurrents + 
        BalanceSplittingCurrents))];
  (*Print["The matrices B and K are: \n",MatrixForm/@
  CoefficientArrays[Kirchhoff,vars=RandomSample@Join[js,jts]],
  "\nThe order of the variables is \n",
  vars];*)(*Outgoing currents at entrances*)
  RuleEntryOut = (jvars[#] -> 0) & /@ (AtTail /@ InEdges);(*Rule*)
  RuleEntryIn = Join[RuleEntryIn, RuleEntryOut];
  (*Not necessary to replace RuleEntryIn in InitRules*)
  AssociateTo[InitRules, RuleEntryIn];
  (*Include Gathering currents information in the rules*){TrueEq, 
    InitRules} = 
   CleanEqualities[{EqBalanceGatheringCurrents /. 
      RuleBalanceGatheringCurrents, InitRules}];
  ExitNeighbors = 
   IncidenceList[AuxiliaryGraph, OutwardVertices /@ ExitVertices];
  Print[ExitNeighbors];
  Print[OutEdges];
  (*T ODO this seems to be defined already:
  OutEdges*)(*Incoming currents at the exits are zero*)
  RuleExitCurrentsIn = ExitCurrents[jvars] /@ ExitNeighbors;(*Rule*)
  Kirchhoff = Kirchhoff /. Join[RuleExitCurrentsIn, RuleEntryOut];
  {BM, KM} = 
   CoefficientArrays[Kirchhoff, 
    vars = Variables[Kirchhoff /. Equal -> Plus]];
  Print["The matrices B and K are: \n", MatrixForm /@ {-BM, KM}, 
   "\nThe order of the variables is \n", vars];
  AssociateTo[InitRules, RuleExitCurrentsIn];
  (*Exit values at exit vertices*)
  RuleExitValues = 
   ExitRules[uvars, ExitCosts] /@ 
    ExitNeighbors;(*Rule*)(*Not necessary to replace RuleExitValues \
in InitRules,there are no us up to now.*)Print[RuleExitValues];
  AssociateTo[InitRules, RuleExitValues];
  (*The value function on the auxiliary edges is constant and equal \
to the exit cost.*)
  EqValueAuxiliaryEdges = 
   And @@ ((uvars[AtHead[#]] == uvars[AtTail[#]]) & /@ 
      Join[InEdges, OutEdges]);(*Equal*)
  RuleValueAuxiliaryEdges = (uvars[AtTail[#]] -> uvars[AtHead[#]]) & /@
     Join[InEdges, OutEdges];(*Equal*)
  Print["D2E: CleanEqualities for the values at the auxiliary \
edges"];
  {TrueEq, InitRules} = 
   CleanEqualities[{EqValueAuxiliaryEdges, InitRules}];
  (*Print[
  RuleValueAuxiliaryEdges];*)(*Infinite switching costs here prevent \
the network from sucking agents from the exits.*)
  OutRules = 
   Rule[#, 
      Infinity] & /@ (Outer[Flatten[{AtTail[#1], #2}] &, OutEdges, 
       EL] // Flatten[#, 1] &);
  InRules = 
   Rule[#, 
      Infinity] & /@ (Outer[{#2[[2]], #1, #2} &, 
       IncidenceList[FG, #] & /@ EntranceVertices, InEdges] // 
      Flatten[#, 2] &);
  AssociateTo[SwitchingCosts, Association[OutRules]];
  AssociateTo[SwitchingCosts, Association[InRules]];
  Print["D2E: Assembled most elements of the system"];
  (*Switching condition equations*)
  EqSwitchingByVertex = 
   Transu[uvars, SwitchingCosts] /@ TransitionsAt[FG, #] & /@ VL;
  (*EqSwitchingByVertex=BooleanConvert[Reduce[Reduce@#,Reals],"CNF"]&/@
  EqSwitchingByVertex;*)(*Maybe this is good with nonzeroswitching \
costs.*)EqSwitchingByVertex = 
   BooleanConvert[Reduce[#, Reals], "CNF"] & /@ EqSwitchingByVertex;
  EqSwitchingByVertex = DeleteCases[EqSwitchingByVertex, True];
  Print["D2E: CleanEqualities for the switching conditions on each \
vertex: \n", Reduce@EqSwitchingByVertex];
  {EqSwitchingConditions, InitRules} = 
   CleanEqualities[{EqSwitchingByVertex, InitRules}];
  EqCompCon = 
   And @@ Compu[jtvars, uvars, SwitchingCosts] /@ AllTransitions;(*Or*)
  Print["D2E: CleanEqualities for the complementary conditions (given \
the rules)"];
  {EqCompCon, InitRules} = CleanEqualities[{EqCompCon, InitRules}];
  (*Default in a,function,
  corresponds to the classic critical congestion case*)
  a (*[SignedCurrents[edge],edge]:*)= 
   Lookup[Data, "a", Function[{j, edge}, j]];
  MinimalTimeRhs = 
   Flatten[-a[SignedCurrents[#], #] + SignedCurrents[#] & /@ BEL];
  Nlhs = 
   Flatten[
    uvars[AtHead[#]] - uvars[AtTail[#]] + SignedCurrents[#] & /@ 
     BEL];
  EqMinimalTime = 
   And @@ (MapThread[(#1 == #2) &, {Nlhs, MinimalTimeRhs}]);
  Nrhs = 
   Flatten[-Cost[SignedCurrents[#], #] + SignedCurrents[#] & /@ 
     BEL];(*one possible cost is IntM*)
  Print["D2E: CleanEqualities for the balance conditions in terms of \
(mostly) transition currents"];
  {TrueEq, InitRules} = 
   CleanEqualities[{EqBalanceSplittingCurrents, InitRules}];
  AllOr = EqCurrentCompCon && EqTransitionCompCon && EqCompCon;
  AllOr = AllOr /. InitRules;
  AllOr = BooleanConvert[Simplify /@ AllOr, "CNF"];
  {AllOr, InitRules} = CleanEqualities[{AllOr, InitRules}];
  EqPosCon = EqPosCon /. InitRules;
  EqSwitchingConditions = EqSwitchingConditions /. InitRules;
  EqSwitchingConditions = Simplify[EqSwitchingConditions];
  AllOr = AllOr && Select[EqSwitchingConditions, Head[#] === Or &];
  AllOr = BooleanConvert[AllOr, "CNF"];
  AllIneq = 
   EqPosCon && Select[EqSwitchingConditions, Head[#] =!= Or &];
  AllIneq = Simplify[AllIneq];
  (*AllOr=BooleanConvert[Simplify[#,AllIneq]&/@AllOr,"CNF"];
  {EqAllAll,InitRules}=CleanEqualities[{AllOr&&AllIneq,InitRules}];*)
  EqCriticalCase = And @@ ((# == 0) & /@ Nlhs);(*Equal*)
  EqNonCritical = 
   And @@ (MapThread[
       Equal[#1, #2] &, {Nlhs, costpluscurrents /. InitRules}]) /. 
    AssociationThread[us, us /. InitRules];
  RuleNonCritical1 = Solve[EqNonCritical, Reals] // Quiet;
  (*The operations below are commutative! Print[RuleNonCritical1/. 
  AssociationThread[costpluscurrents,0&/@costpluscurrents]];
  Print[Solve[EqNonCritical/. AssociationThread[costpluscurrents,0&/@
  costpluscurrents],Reals]//Quiet];*)(*Print[
  EqCriticalCase/.AssociationThread[us,us/.InitRules]];
  Print[EqNonCritical/. AssociationThread[costpluscurrents,0&/@
  costpluscurrents]];
  Print["right?\n",(EqCriticalCase/.AssociationThread[us,
  us/.InitRules])/.RuleNonCritical1];
  Print["right?\n",(EqCriticalCase/. AssociationThread[
  costpluscurrents,0&/@costpluscurrents])/.(RuleNonCritical1/. 
  AssociationThread[costpluscurrents,0&/@costpluscurrents])];
  Print[EqNonCritical/. RuleNonCritical1];*)
  Print["D2E: CleanEqualities for the (non) critical case \
equations"];
  {TrueEq, RuleNonCritical} = 
   CleanEqualities[{EqNonCritical, InitRules}];
  (*Print[Expand/@RuleNonCritical];*)
  EqCriticalCase = 
   EqNonCritical /. 
    AssociationThread[costpluscurrents, 0 & /@ costpluscurrents];
  (*Print[EqCriticalCase];*)
  RulesCriticalCase1 = 
   Expand /@ (RuleNonCritical1 /. 
      AssociationThread[costpluscurrents, 0 & /@ costpluscurrents]);
  RulesCriticalCase = 
   Expand /@ (RuleNonCritical /. 
      AssociationThread[costpluscurrents, 0 & /@ costpluscurrents]);
  (*Here we change the pourpose of costpluscurrents*)
  costpluscurrents = AssociationThread[costpluscurrents, Nrhs];
  (*Print[Solve[EqCriticalCase,js]];
  (*RulesCriticalCaseJs=Association@First[Solve[
  EqCriticalCase/.RuleEntryIn,js]];*)Print[
  "CleanEqualities for the critical case equations"];
  {TrueEq,RulesCriticalCase}=CleanEqualities[{EqCriticalCase,
  InitRules}];
  Print["Expanding critical rules..."];
  RulesCriticalCase=Expand/@RulesCriticalCase;
  Print[RulesCriticalCase];*)Print["D2E: D2E is finished!"];
  Print[$Context];
  Print[Names[$Context <> "*"]];
  ttt = 1;
  Print[Names[$Context <> "*"]];
  Join[Data, 
   Association[(*Graph structure*)"BG" -> BG, "InEdges" -> InEdges, 
    "OutEdges" -> OutEdges, "FG" -> FG,(*variables*)"jvars" -> jvars, 
    "jtvars" -> jtvars, "uvars" -> uvars, 
    "costpluscurrents" -> costpluscurrents, 
    "jays" -> SignedCurrents,(*equations*)(*complementarity*)
    "AllOr" -> AllOr,(*union of all complementarity conditions*)
    "AllIneq" -> AllIneq, "EqPosJs" -> EqPosJs, 
    "EqPosJts" -> EqPosJts, "EqPosCon" -> EqPosCon, 
    "EqCurrentCompCon" -> EqCurrentCompCon, 
    "EqTransitionCompCon" -> 
     EqTransitionCompCon,(*linear equations (and inequalities)*)
    "EqAllAll" -> EqAllAll, "BoundaryRules" -> InitRules, 
    "InitRules" -> InitRules, 
    "RulesCriticalCase" -> RulesCriticalCase, 
    "RulesCriticalCase1" -> 
     RulesCriticalCase1,(*"RulesCriticalCaseJs"->RulesCriticalCaseJs,*)
    "RuleExitValues" -> RuleExitValues, "RuleEntryIn" -> RuleEntryIn, 
    "Nlhs" -> Nlhs, "MinimalTimeRhs" -> MinimalTimeRhs, 
    "EqCriticalCase" -> EqCriticalCase, 
    "EqNonCritical" -> EqNonCritical,(*"EqNonCritical1"->
    EqNonCritical1,*)"RuleNonCritical" -> RuleNonCritical, 
    "RuleNonCritical1" -> RuleNonCritical1, 
    "EqSwitchingConditions" -> And @@ EqSwitchingByVertex, 
    "EqValueAuxiliaryEdges" -> EqValueAuxiliaryEdges, 
    "EqBalanceSplittingCurrents" -> EqBalanceSplittingCurrents, 
    "EqBalanceGatheringCurrents" -> EqBalanceGatheringCurrents, 
    "Nrhs" -> Nrhs, "B" -> BM, "K" -> KM, "NoDeadEnds" -> NoDeadEnds, 
    "NoDeadStarts" -> NoDeadStarts, "EqEntryIn" -> EqEntryIn, 
    "vars" -> vars]]]


(*"EntranceVertices"->EntranceVertices,*)
\
(*"InwardVertices"->InwardVertices,*)
(*"ExitVertices"->ExitVertices,*)
\
(*"OutwardVertices"->OutwardVertices,*)
\
(*"AuxiliaryGraph"->AuxiliaryGraph,*)
\
(*"VL"->VL,"EL"->EL,"BEL"->BEL,"FVL"->FVL,*)
\
(*"AllTransitions"->AllTransitions,"NoDeadEnds"->NoDeadEnds,\
"NoDeadStarts"->NoDeadStarts,*)
(*"jargs"->jargs,"js"->js,*)
\
(*"jts"->jts,*)
(*"uargs"->uargs,"us"->us,*)
\
(*"SwitchingCosts"->SwitchingCosts,*)
(*"OutRules"->OutRules,*)
\
(*"InRules"->InRules,*)
(*"EntryArgs"->EntryArgs,*)
\
(*"EntryDataAssociation"->EntryDataAssociation,*)
\
(*"ExitCosts"->ExitCosts,*)
(*"EqCompCon"->EqCompCon,*)
\
(*"EqPosCon"->EqPosCon,"EqEntryIn"->EqEntryIn,"EqExitValues"->\
EqExitValues,"EqValueAuxiliaryEdges"->EqValueAuxiliaryEdges,"AllIneq"->\
AllIneq,*)



End[]