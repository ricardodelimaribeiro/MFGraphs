(*Wolfram Language package*)
ConsistentSwithingCosts::usage = "ConsistentSwithingCosts[switching \
costs][one switching cost]  returns the condition for this switching \
cost to satisfy the triangle inequality";

ConsistentSwithingCosts[sc_][{a_, b_, c_} -> S_] := 
  Module[{origin, destination, bounds}, 
   origin = Cases[sc, HoldPattern[{a, b, _} -> _]];
   origin = DeleteCases[origin, {a, b, c} -> S];
   If[origin =!= {},
    destination = {a, Part[First[#], 3], c} & /@ origin;
    bounds = ((S <= 
          Last[#] + Association[sc][{a, Part[First[#], 3], c}]) & /@ 
       origin);
    And @@ bounds, True]];

IsSwitchingCostConsistent::usage = "IsSwitchingCostConsistent[List of \
switching costs] is True if all switching costs satisfy the triangle \
inequality"
IsSwitchingCostConsistent[SC_] := 
 And @@ Simplify[ConsistentSwithingCosts[SC] /@ SC]

AtHead[DirectedEdge[a_, b_]] := {b, DirectedEdge[a, b]}

AtTail[DirectedEdge[a_, b_]] := {a, DirectedEdge[a, b]}

TransitionsAt[G_, k_] := 
 Prepend[#, k] & /@ Permutations[IncidenceList[G, k], {2}]

RoundValues[x_?NumberQ] := Round[x, 10^-10]

RoundValues[Rule[a_, b_]] := Rule[a, RoundValues[b]]

RoundValues[x_List] := RoundValues /@ x

RoundValues[x_Association] := RoundValues /@ x

triple2path[{a_, b_, c_}, G_] := 
 Module[{EL = EdgeList[G]}, 
  If[SubsetQ[AdjacencyList[G, b], {a, c}], 
   If[MemberQ[EL, DirectedEdge[a, b]], 
    If[MemberQ[EL, DirectedEdge[b, c]], {b, DirectedEdge[a, b], 
      DirectedEdge[b, c]}, {b, DirectedEdge[a, b], 
      DirectedEdge[c, b]}], 
    If[MemberQ[EL, DirectedEdge[b, a]], 
     If[MemberQ[EL, DirectedEdge[b, c]], {b, DirectedEdge[b, a], 
       DirectedEdge[b, c]}, {b, DirectedEdge[b, a], 
       DirectedEdge[c, b]}]]], 
   StringJoin["\nThere is no path from ", ToString[a], " to ", 
    ToString[b], " to ", ToString[c]]]]

path2triple[{a_, b_ \[DirectedEdge] c_, d_ \[DirectedEdge] e_} -> 
    S_] := {Select[{b, c}, # =!= a &], a, Select[{d, e}, # =!= a &], 
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

IncomingEdges[FG_][k_] := {k, #1} & /@ IncidenceList[FG, k];

OtherWay[{c_, DirectedEdge[a_, b_]}] := {If[c === a, b, a], 
  DirectedEdge[a, b]}

OutgoingEdges[FG_][k_] := 
  OtherWay /@ ({k, #} & /@ IncidenceList[FG, k]);

ExitRules[uvars_, ExitCosts_][a_ \[DirectedEdge] b_] := 
  Total[uvars /@ {{b, DirectedEdge[a, b]}}] -> ExitCosts[b];

ExitCurrents[jvars_][a_ \[DirectedEdge] b_] := 
  jvars@{a, DirectedEdge[a, b]} -> 0;

Transu[uvars_, SwitchingCosts_][{v_, edge1_, edge2_}] := 
  uvars[{v, edge1}] <= 
   uvars[{v, edge2}] + SwitchingCosts[{v, edge1, edge2}];

Compu[jtvars_, uvars_, SwitchingCosts_][{v_, edge1_, 
    edge2_}] := (jtvars[{v, edge1, edge2}] == 0) || 
   uvars[{v, edge2}] - uvars[{v, edge1}] + 
     SwitchingCosts[{v, edge1, edge2}] == 0;

Data2Equations::usage = "Data2Equations[Data] returns the equations, \
inequalities, and alternatives associated to the Data"
Data2Equations[Data_Association] := 
  Module[{VL, AM, EVC, EVTC, SC, BG, EntranceVertices, InwardVertices,
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
    consistentCosts, costpluscurrents, EqGeneral}, 
   VL = Lookup[Data, "Vertices List", {}];
   AM = Lookup[Data, "Adjacency Matrix", {}];
   EVC = Lookup[Data, "Entrance Vertices and Currents", {}];
   EVTC = Lookup[Data, "Exit Vertices and Terminal Costs", {}];
   SC = Lookup[Data, "Switching Costs", {}];
   (****Graph stuff****)
   BG = AdjacencyGraph[VL, AM, VertexLabels -> "Name", 
     DirectedEdges -> True];
   EntranceVertices = First /@ EVC;
   ExitVertices = First /@ EVTC;
   Clear["en*", "ex*"];
   (*InwardVertices defines auxiliary vertices for the entrance \
vertices*)
   InwardVertices = 
    AssociationMap[Symbol["en" <> ToString[#]] &, EntranceVertices];
   OutwardVertices = 
    AssociationMap[Symbol["ex" <> ToString[#]] &, ExitVertices];
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
   EL = EdgeList[FG];
   BEL = EdgeList[BG];
   FVL = VertexList[FG];
   (*arguments*)
   jargs = Flatten[#, 1] &@({AtTail@#, AtHead@#} & /@ EL);
   uargs = jargs;
   (*AllTransitions=TransitionsAt[BG,#]&/@VL//
   Catenate(*at vertex from first edge to second edge*);*)
   AllTransitions = 
    TransitionsAt[FG, #] & /@ FVL // 
     Catenate(*at vertex from first edge to second edge*);
   EntryArgs = 
    AtHead /@ ((EdgeList[AuxiliaryGraph, 
           DirectedEdge[_, #]] & /@ (First /@ EVC)) // 
       Flatten[#, 1] &);
   ZeroRun = 
    AssociationMap[(0 &) &][AtHead /@ Join[InEdges, OutEdges]];
   EntryDataAssociation = 
    RoundValues@AssociationThread[EntryArgs, Last /@ EVC];
   ExitCosts = 
    AssociationThread[OutwardVertices /@ (First /@ EVTC), 
     Last /@ EVTC];
   (*variables*)
   js = Table[Symbol["j" <> ToString[k]], {k, 1, Length@jargs}];
   jvars = AssociationThread[jargs, js];
   jts = 
    Table[
     Symbol["jt" <> ToString[k]], {k, 1, Length@AllTransitions}];
   jtvars = AssociationThread[AllTransitions, jts];
   us = Table[Symbol["u" <> ToString[k]], {k, 1, Length@uargs}];
   uvars = AssociationThread[uargs, us];
   SignedCurrents = 
    AssociationMap[jvars[AtHead[#]] - jvars[AtTail[#]] &, BEL];
   SwitchingCosts = AssociationMap[0 &, AllTransitions];
   AssociateTo[SwitchingCosts, 
    AssociationThread[triple2path[Take[#, 3], FG] & /@ SC, 
     Last[#] & /@ SC]];
   LargeCases = 
    Join[{Last[#], __, #} & /@ InEdges, {First[#], #, __} & /@ 
      OutEdges];
   LargeSwitchingTransitions = 
    Cases[AllTransitions, #] & /@ LargeCases // Flatten[#, 1] &;
   AssociateTo[SwitchingCosts, 
    AssociationMap[1000000 &, LargeSwitchingTransitions]];
   consistentCosts = 
    IsSwitchingCostConsistent[Normal@SwitchingCosts];
   Which[consistentCosts === False, 
   	Print["Switching costs are inconsistent!"];
    Return[ $Failed, Module], 
    consistentCosts =!= True, 
    Print["Switching costs conditions are ", consistentCosts](*, True, 
    Print["Switching costs are consistent"]*)];
   CostArgs = 
    Join[AssociationMap[Identity &, Normal@jargs], ZeroRun, 
     AssociationMap[(10^(-4) &) &, Normal@Keys@SwitchingCosts], 
     AssociationMap[(10^4 &) &, LargeSwitchingTransitions]];
   EqPosJs = And @@ (# >= 0 & /@ Join[jvars]);(*Inequality*)
   EqPosJts = And @@ (# >= 0 & /@ Join[jtvars]);(*Inequality*)
   EqCurrentCompCon = And @@ (CurrentCompCon[jvars] /@ EL);(*Or*)
   EqTransitionCompCon = 
    And @@ ((Sort /@ TransitionCompCon[jtvars] /@ AllTransitions) // 
       Union);(*Or*)(*Balance Splitting Currents in the full graph*)
   NoDeadEnds = IncomingEdges[FG] /@ VL // Flatten[#, 1] &;
   BalanceSplittingCurrents = ((jvars[#] - 
         Total[jtvars /@ CurrentSplitting[AllTransitions][#]]) & /@ 
      NoDeadEnds);
   EqBalanceSplittingCurrents = 
    Simplify /@ (And @@ ((# == 0) & /@ 
         BalanceSplittingCurrents));(*Equal*)(*Gathering currents in \
the inside of the basic graph*)
   NoDeadStarts = OutgoingEdges[FG] /@ VL // Flatten[#, 1] &;
   RuleBalanceGatheringCurrents = 
    Association[(jvars[#] -> 
         Total[jtvars /@ CurrentGathering[AllTransitions][#]]) & /@ 
      NoDeadStarts];(*Rule*)(*get equations for the exit currents at \
the entry vertices*)
   BalanceGatheringCurrents = ((-jvars[#] + 
         Total[jtvars /@ CurrentGathering[AllTransitions][#]]) & /@ 
      NoDeadStarts);
   EqBalanceGatheringCurrents = 
    Simplify /@ (And @@ (# == 0 & /@ BalanceGatheringCurrents));
   (*Incoming currents*)
   EqEntryIn = (jvars[#] == EntryDataAssociation[#]) & /@ (AtHead /@ 
       InEdges);(*List of Equals*)(*Outgoing currents at entrances*)
   RuleEntryOut = 
    Association[(jvars[#] -> 0) & /@ (AtTail /@ InEdges)];(*Rule*)
   RuleExitCurrentsIn = 
    Association[
     ExitCurrents[jvars] /@ 
      OutEdges];(*Rule*)(*Exit values at exit vertices*)
   RuleExitValues = 
    Association[
     ExitRules[uvars, ExitCosts] /@ 
      OutEdges];(*Rule*)(*The value function on the auxiliary edges \
is constant and equal to the exit cost.*)
   EqValueAuxiliaryEdges = 
    And @@ ((uvars[AtTail[#]] == uvars[AtHead[#]]) & /@ 
       Join[InEdges, 
        OutEdges]);(*Equal*)(*use ToRules to get the rules*)
   OutRules = 
    Rule[#, 
       Infinity] & /@ (Outer[Flatten[{AtTail[#1], #2}] &, OutEdges, 
        EL] // Flatten[#, 1] &);
   InRules = 
    Rule[#, 
       Infinity] & /@ (Outer[{#2[[2]], #1, #2} &, 
        IncidenceList[FG, #] & /@ EntranceVertices, InEdges] // 
       Flatten[#, 2] &);
   EqSwitchingByVertex = 
    Transu[uvars, SwitchingCosts] /@ TransitionsAt[FG, #] & /@ VL;
   EqSwitchingByVertex = (And @@ #) & /@ EqSwitchingByVertex;
   EqCompCon = 
    And @@ 
     Compu[jtvars, uvars, SwitchingCosts] /@ AllTransitions;(*Or*)
   Nlhs = 
    Flatten[
     uvars[AtHead[#]] - uvars[AtTail[#]] + SignedCurrents[#] & /@ BEL];
   Nrhs = 
    Flatten[SignedCurrents[#] - IntM[SignedCurrents[#], #] & /@ BEL];
  (*stuff to solve the general case faster*)
   costpluscurrents = 
    Table[Symbol["cpc" <> ToString[k]], {k, 1, Length@BEL}];
    (*Nrhs = costpluscurrents;*)
   EqGeneral = 
    And @@ (MapThread[Equal, {Nlhs, costpluscurrents}]);
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
     CostArgs, Nrhs, costpluscurrents, EqGeneral};
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
     "Nrhs", "costpluscurrents", "EqGeneral"};
   Join[Data, AssociationThread[ModuleVarsNames, ModuleVars]]];

NumberVectorQ[j_] := And @@ (NumberQ /@ j);

GetKirchhoffMatrix::usage = "GetKirchhoffMatrix[d2e] returns the \
Kirchhoff matrix, entry current vector, (critical congestion) cost \
function, and the variables order."
GetKirchhoffMatrix[Eqs_] := 
  Module[{Kirchhoff, EqEntryIn = Lookup[Eqs, "EqEntryIn", True], 
    BalanceGatheringCurrents = 
     Lookup[Eqs, "BalanceGatheringCurrents", {}], 
    BalanceSplittingCurrents = 
     Lookup[Eqs, "BalanceSplittingCurrents", {}], 
    RuleExitCurrentsIn = Lookup[Eqs, "RuleExitCurrentsIn", {}], 
    RuleEntryOut = Lookup[Eqs, "RuleEntryOut", {}], BM, KM, vars, 
    CostArgs = Lookup[Eqs, "CostArgs", <||>], 
    jvars = Lookup[Eqs, "jvars", {}], 
    jtvars = Lookup[Eqs, "jtvars", {}], cost}, 
   Kirchhoff = 
    Join[
     EqEntryIn, (# == 0 & /@ (BalanceGatheringCurrents + 
         BalanceSplittingCurrents))];
   Kirchhoff = Kirchhoff /. Join[RuleExitCurrentsIn, RuleEntryOut];
   {BM, KM} = 
    CoefficientArrays[Kirchhoff, 
     vars = Variables[Kirchhoff /. Equal -> Plus]];
   cost = 
    Function[currents, 
     MapThread[#1[#2] &, {KeyMap[Join[jvars, jtvars]][CostArgs] /@ 
        vars, currents}]];
   {-BM, KM, cost, vars}];

MFGPreprocessing[Eqs_] := 
  Module[{InitRules, RuleBalanceGatheringCurrents, EqEntryIn, 
    RuleEntryOut, RuleExitCurrentsIn, RuleExitValues, 
    EqValueAuxiliaryEdges, EqSwitchingByVertex, EqCompCon, 
    EqBalanceSplittingCurrents, EqCurrentCompCon, EqTransitionCompCon,
     EqPosJs, EqPosJts, ModuleVarsNames, ModulesVars, NewSystem, 
    Rules, EqGeneral}, 
   RuleBalanceGatheringCurrents = 
    Lookup[Eqs, "RuleBalanceGatheringCurrents", $Failed];
   EqEntryIn = Lookup[Eqs, "EqEntryIn", $Failed];
   RuleEntryOut = Lookup[Eqs, "RuleEntryOut", $Failed];
   RuleExitCurrentsIn = Lookup[Eqs, "RuleExitCurrentsIn", $Failed];
   RuleExitValues = Lookup[Eqs, "RuleExitValues", $Failed];
   EqCompCon = Lookup[Eqs, "EqCompCon", $Failed];
   EqCurrentCompCon = Lookup[Eqs, "EqCurrentCompCon", $Failed];
   EqTransitionCompCon = Lookup[Eqs, "EqTransitionCompCon", $Failed];
   EqPosJs = Lookup[Eqs, "EqPosJs", $Failed];
   EqPosJts = Lookup[Eqs, "EqPosJts", $Failed];
   EqGeneral =Lookup[Eqs, "EqGeneral", $Failed];
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
   EqSwitchingByVertex = 
    And @@ (Simplify /@ ( EqSwitchingByVertex /. InitRules));
   EqBalanceSplittingCurrents = EqBalanceSplittingCurrents /. InitRules;
   EqValueAuxiliaryEdges = EqValueAuxiliaryEdges /. InitRules;
   Rules = 
    First@Solve[EqBalanceSplittingCurrents && EqValueAuxiliaryEdges] //
      Quiet;
   InitRules = Join[InitRules /. Rules, Association@Rules];
   NewSystem = 
    MapThread[
     And, {{True, EqPosJts && EqPosJs, 
       EqCurrentCompCon && EqTransitionCompCon && EqCompCon}, 
      Sys2Triple[EqSwitchingByVertex]}];
   {NewSystem, InitRules} = TripleClean[{NewSystem, InitRules}];
   NewSystem[[1]] = EqGeneral;
   {NewSystem, InitRules} = TripleClean[{NewSystem, InitRules}];
   (*sometimes Reduce does not reduce, but if we simplify each term it does.*)
   NewSystem[[3]] = Simplify /@ NewSystem[[3]];
   ModuleVarsNames = {"InitRules", "NewSystem"};
   ModulesVars = {InitRules, NewSystem} ;
   Join[Eqs, AssociationThread[ModuleVarsNames, ModulesVars]]
   ];

CriticalCongestionSolver::usage = 
  "CriticalCongestionSolver[Eqs] returns an association with rules to \
the solution to the critical congestion case";
CriticalCongestionSolver[$Failed] := $Failed

CriticalCongestionSolver[Eqs_] := 
  Module[{PreEqs, js, AssoCritical}, 
   PreEqs = MFGPreprocessing[Eqs];
   js = Lookup[PreEqs, "js",$Failed];
   AssoCritical = MFGSystemSolver[PreEqs][AssociationThread[js, 0 js]];
   Join[PreEqs, Association["AssoCritical"-> AssoCritical]]
   ];

MFGSystemSolver::usage = "MFGSystemSolver[Eqs][edgeEquations] returns \
an association with rules to the solution";
MFGSystemSolver[Eqs_][approxJs_] := 
  Module[{NewSystem, InitRules, pickOne, vars, System, Ncpc,
  	costpluscurrents}, 
   InitRules = Lookup[Eqs, "InitRules", $Failed];
   NewSystem = Lookup[Eqs, "NewSystem", $Failed];
   costpluscurrents = Lookup[Eqs, "costpluscurrents", $Failed];
   Ncpc = RoundValues[Expand/@(costpluscurrents /. approxJs)];
   InitRules = Expand/@(InitRules /. Ncpc);
   NewSystem = NewSystem /. Ncpc;
   {NewSystem, InitRules} = FinalClean[{NewSystem, InitRules}];
   System = And @@ NewSystem;
   Which[System === False, Print["MFGSS: There is no solution"], 
    System =!= True, 
    NewSystem = Reduce[System, Reals];
    {NewSystem, InitRules} = FinalClean[{Sys2Triple[NewSystem], InitRules}];
    NewSystem = Reduce[And @@ NewSystem, Reals];
    (*not checking if NewSystem is not True...*)
    Print["MFGSS: Multiple solutions: ", NewSystem];
    vars = 
     Select[Join[Eqs["us"], Eqs["js"], Eqs["jts"]], 
      Not[FreeQ[NewSystem, #]] &];
    (*Have to pick one so that all the currents have numerical values*)
    pickOne = 
     Association@
      First@
       FindInstance[NewSystem && And @@ ((# >(*=*) 0 )& /@ vars), vars, 
        Reals];
    Print["\tPicked one value for the variable(s) ", vars];
    InitRules = Expand /@ Join[InitRules /. pickOne, pickOne]
    ];
    InitRules
    ];

FinalStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, 
   rules_}] := {{EE, NN, OR}, rules}

FinalStep[{{EE_, NN_, OO_}, rules_}] := 
  Module[{NewSystem, newrules, sorted, time}, 
  	{NewSystem, newrules} = TripleClean[{{EE, NN, OO}, rules}];
   If[Part[NewSystem, 3] === True, sorted = True, 
    sorted = ReverseSortBy[Part[NewSystem, 3], Simplify`SimplifyCount]];
    {time, NewSystem} = 
    AbsoluteTiming[ZAnd[And @@ Take[NewSystem, {1, 2}], 
       If[Part[NewSystem, 3] === True, True, sorted]]];
   Print["Iterative boolean convert took ", time, " seconds to terminate"];
   If[Head[NewSystem] === Or, NewSystem = Sort /@ NewSystem];
   NewSystem = Sys2Triple[NewSystem];
   {NewSystem, newrules}];

FinalClean[{{EE_, NN_, OR_}, rules_}] := ((*Print["finalclean: ",{{EE,
  NN,OR},rules}];*)FixedPoint[FinalStep, {{EE, NN, OR}, rules}])

Sys2Triple[True] = Table[True, 3]

Sys2Triple[False] = Table[False, 3]

Sys2Triple[system_] := 
  Which[Head[system] === And, 
   Module[{groups, EE, OR, NN}, 
    groups = GroupBy[List @@ system, Head[#] === Equal &];
    EE = And @@ Lookup[groups, True, {}];
    groups = GroupBy[Lookup[groups, False, {}], Head[#] === Or &];
    OR = And @@ Lookup[groups, True, True];
    NN = And @@ Lookup[groups, False, True];
    {EE, NN, OR}], Head[system] === Equal, {system, True, True}, 
   Head[system] === Or, {True, True, system}, 
   True, {True, system, True}];

TripleStep[{{EEs_, NNs_, ORs_}, rules_List}] := 
 TripleStep[{{EEs, NNs, ORs}, Association@rules}]

TripleStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, 
   rules_}] := {{EE, NN, OR}, rules}

TripleStep[{{EE_, NN_?TrueQ, OR_?TrueQ}, rules_Association}] := 
  Module[{newrules = {}}, newrules = First@Solve[EE] // Quiet;
   newrules = Join[rules /. newrules, Association@newrules];
   {{True, NN, OR}, Expand /@ newrules}];

TripleStep[{{EEs_, NNs_, ORs_}, rules_Association}] := 
  Module[{EE = EEs /. rules, NN = Simplify /@ (NNs /. rules), 
    OR = Simplify /@ (ORs /. rules), NNE, NNO, ORE, ORN, bool, 
    newrules = {}}, bool = EE && NN && OR;
   NN = Simplify[NN];
   {NNE, NN, NNO} = Sys2Triple[NN];
   {ORE, ORN, OR} = Sys2Triple[OR];
   EE = EE && NNE && ORE;
   If[EE =!= True && EE =!= False, 
    newrules = First@Solve[EE] // Quiet];
   newrules = 
    Expand /@ Join[rules /. newrules, Association@newrules];
   (*in=Simplify/@((And[EEs,NNs,ORs]/.rules)&&(And@@Equal@@@Normal@
   rules));
   out=Simplify/@((And[EE,NN,OR]/.newrules)&&(And@@Equal@@@Normal@
   newrules));
   Print["Checking..."];
   Print[in];
   Print[out];*)(*in=Simplify/@(And[EEs,NNs,ORs]/.newrules);
   out=Simplify/@(And[EE,NN,OR]/.newrules);
   Print[in];
   Print[out];
   Print["Checking (replacing updated rules in both systems)..."];(**)
   If[Reduce[Implies[in,out],Reals]===True,Print["Ok!"],Print[
   "Not Ok..."]];*){{EE, NN, OR}, Expand /@ newrules}];


TripleClean[{{EE_, NN_, OR_}, rules_}] := 
 FixedPoint[TripleStep, {{EE, NN, OR}, rules}]

NewReduce[s_Or] := Reduce[s, Reals]

NewReduce[x_] := (Print["NewReduce[", x, "] = ", x];
  x)

(*TODO NewReduce:is there a problem here?*)

NewReduce[system_And] := 
 Module[{result, 
   groups = 
    GroupBy[List @@ system, Head[#] === Or || Head[#] === Equal &], 
   subgroups, EE, NN, OO, sorted}, 
  subgroups = GroupBy[groups[True], Head[#] === Equal &];
  NN = Lookup[groups, False, True];
  EE = Lookup[subgroups, True, True];
  OO = Lookup[subgroups, False, True];
  If[OO === True, result = ZAnd[And @@ NN, And @@ EE], 
   sorted = SortBy[And @@ OO, Simplify`SimplifyCount];
   result = ZAnd[And @@ NN, (And @@ EE) && sorted]];
  If[result =!= False, result = result // DeleteDuplicates];
  result]
ZAnd::usage =
"ZAnd[processed, unprocessed] does something...";

ZAnd[_, False] := False

ZAnd[False, _] := False

ZAnd[xp_, True] := xp

ZAnd[xp_, eq_Equal] := 
 With[{sol = Solve[eq]}, 
  If[sol === {}, 
   False, (xp /. First@sol) && And @@ (First@sol /. Rule -> Equal) // 
    Simplify]]

ZAnd[xp_, andxp_And] := 
 With[{fst = First[andxp], rst = Rest[andxp]}, 
  Which[Head[fst] === Or, ReZAnd[xp, rst] /@ fst // RemoveDuplicates, 
   True, ReZAnd[xp, rst, fst]]]

ZAnd[xp_, orxp_Or] := ((*Print["Do we need removeduplicates?\n",(ZAnd[
  xp,#]&/@orxp)];*)(ZAnd[xp, #] & /@ orxp) // RemoveDuplicates)

ZAnd[xp_, leq_] := Simplify[xp && leq]

(*ZAnd[xp_,geq_GreaterEqual]:=Simplify[xp&&geq] \
ZAnd[xp_,ineq_Inequality]:=Simplify[xp&&ineq]*)


(*Operator form of ReZAnd*)
ReZAnd[xp_, rst_] := ReZAnd[xp, rst, #] &

ReZAnd[xp_, rst_, fst_Equal] := 
 Module[{fsol = First@Solve@fst // Quiet, newrst, newxp}, 
  newrst = ReplaceSolution[rst, fsol];
  newxp = xp /. fsol;
  ZAnd[newxp && fst, newrst]]

(*ReZAnd[xp_,rst_,fst_And]:=ReZAnd[xp,rst]/@fst (*TODO never used???*)*)
\

ReZAnd[xp_, rst_, fst_] := ZAnd[xp && fst, rst]

ReplaceSolution[rst_?BooleanQ, sol_] := rst

ReplaceSolution[rst_, sol_] := Module[{newrst}, newrst = rst /. sol;
  If[Head[newrst] === And, Reduce[#, Reals] & /@ newrst, 
   Reduce[newrst, Reals]]]

(*ReplaceSolution[rst_,sol_]:=Simplify[rst/. sol]*)

RemoveDuplicates[xp_And] := DeleteDuplicates[Sort[xp]];
RemoveDuplicates[xp_Or] := DeleteDuplicates[Sort[xp]];
RemoveDuplicates[xp_] := xp

(****************************************************************)
(****************************************************************)
\
(****************************************************************)
(****************************************************************)
\
(****************************************************************)
(****************************************************************)
\
(****************************************************************)
(****************************************************************)
\
(****************************************************************)