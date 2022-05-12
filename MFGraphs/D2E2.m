(*Wolfram Language package*)
ConsistentSwithingCosts::usage = "ConsistentSwithingCosts[switching \
costs][one switching cost]  returns the condition for this switching \
cost to satisfy the triangle inequality"
ConsistentSwithingCosts[sc_][{a_, b_, c_, S_}] := 
  Module[{origin, destination, bounds}, 
   origin = Cases[sc, {a, b, _, _}];
   destination = Cases[sc, {_, b, c, _}];
   bounds = 
    Outer[Plus, Last /@ origin, Last /@ destination] // Flatten;
   Print[bounds];
   And @@ (# >= S &) /@ bounds];

IsSwitchingCostConsistent::usage = "IsSwitchingCostConsistent[list of \
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
    EqCompCon, Nlhs, ModuleVars, ModuleVarsNames,
    LargeCases, LargeSwitchingTransitions, ZeroRun,
    CostArgs, Nrhs}, 
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
   AssociationMap[ 
     Symbol["en" <> ToString[#]] &, EntranceVertices];
   OutwardVertices = 
    AssociationMap[ 
     Symbol["ex" <> ToString[#]] &, ExitVertices];
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
   (*arguments*)jargs = Flatten[#, 1]& @ ({AtTail@#, AtHead@#} & /@ EL);
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
  (*     Print["inedges: ",InEdges];
   Print["entry: ", EntryArgs];
   Print[AtHead/@InEdges, AtTail/@InEdges];*)
   ZeroRun = AssociationMap[(0&) &][AtHead/@Join[InEdges, OutEdges]];
   (*TODO use largerun to set large costs on the ficticious *)
   (*Print["zerorun: ", ZeroRun];*)
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
    AssociationMap[
     jvars[AtHead[#]] - jvars[AtTail[#]] &, BEL];
   
   (*Print[
   "D2E: Variables are all set"];*)(*Elements of the \
system*)(*Swithing cost is initialized with 0. AssociationThread \
associates the last association!*)(*SwitchingCosts=AssociationThread[
   Join[AllTransitions,triple2path[Take[#,3],BG]&/@SC],Join[0&/@
   AllTransitions,Last[#]&/@SC]];*)
   SwitchingCosts = AssociationMap[0&, AllTransitions];
   
   AssociateTo[SwitchingCosts, AssociationThread[
     triple2path[Take[#, 3], FG] & /@ SC, 
     Last[#] & /@ SC]];
     
   LargeCases = Join[{Last[#], __, #} & /@ InEdges, 
   	{First[#], #, __} & /@ OutEdges];
   LargeSwitchingTransitions = Cases[AllTransitions, #] & /@ LargeCases // 
 Flatten[#, 1] &;

   CostArgs = Join[
   	AssociationMap[Identity & , Normal@jargs], ZeroRun, 
	AssociationMap[(10^(-2)&)&, Normal@Keys@SwitchingCosts],
	AssociationMap[(10^2&) &, 
	LargeSwitchingTransitions]];
   EqPosJs = And @@ (# >= 0 & /@ Join[jvars]);(*Inequality*)
   EqPosJts = And @@ (# >= 0 & /@ Join[jtvars]);(*Inequality*)
   EqCurrentCompCon = And @@ (CurrentCompCon[jvars] /@ EL);(*Or*)
   EqTransitionCompCon = 
    And @@ ((Sort /@ TransitionCompCon[jtvars] /@ AllTransitions) // 
       Union);(*Or*)(*Balance Splitting Currents in the full \
graph*)(*NoDeadEnds=IncomingEdges[BG]/@VL//Flatten[#,1]&;*)
   NoDeadEnds = IncomingEdges[FG] /@ VL // Flatten[#, 1] &;
   BalanceSplittingCurrents = ((jvars[#] - 
         Total[jtvars /@ CurrentSplitting[AllTransitions][#]]) & /@ 
      NoDeadEnds);
   EqBalanceSplittingCurrents = 
    Simplify /@ (And @@ ((# == 0) & /@ 
         BalanceSplittingCurrents));(*Equal*)(*Gathering currents in \
the inside of the basic graph*)(*NoDeadStarts=OutgoingEdges[BG]/@VL//
   Flatten[#,1]&;*)
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
   (*Switching condition equations*)(*EqSwitchingByVertex=Transu[
   uvars,SwitchingCosts]/@TransitionsAt[BG,#]&/@VL;*)
  (* Print["b4: ",SwitchingCosts];
   SwitchingCosts =
    SwitchingCosts /. 0 -> 10^(-10);
   Print["aftr: ",SwitchingCosts];*)
   
   EqSwitchingByVertex = 
    Transu[uvars, SwitchingCosts] /@ TransitionsAt[FG, #] & /@ VL;
    (*Print[EqSwitchingByVertex];*)
   EqSwitchingByVertex = (And @@ #) & /@ EqSwitchingByVertex;
   EqCompCon = 
    And @@ 
     Compu[jtvars, uvars, SwitchingCosts] /@ AllTransitions;(*Or*)
   Nlhs = 
    Flatten[
     uvars[AtHead[#]] - uvars[AtTail[#]] + SignedCurrents[#] & /@ BEL];
   Nrhs = Flatten[SignedCurrents[#]-IntM[SignedCurrents[#],#]& /@ BEL];
   ModuleVars =(*list of all module variables,
    except for ModuleVars*){VL, AM, EVC, EVTC, SC, BG, 
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
     OutRules, InRules, EqSwitchingByVertex, EqCompCon, Nlhs, CostArgs, Nrhs};
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
      "EqSwitchingByVertex", "EqCompCon", "Nlhs", "CostArgs", "Nrhs"};
   Join[Data, AssociationThread[ModuleVarsNames, ModuleVars]]];
GetKirchhoffMatrix::usage =
"GetKirchhoffMatrix[d2e] returns the Kirchhoff matrix, entry current vector, (critical congestion) cost function, and the variables order."
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
    jtvars = Lookup[Eqs, "jtvars", {}]
    }, 
    NumberVectorQ[j_] := And @@ (NumberQ /@ j);
   Kirchhoff = 
    Join[
     EqEntryIn, (# == 0 & /@ (BalanceGatheringCurrents + 
         BalanceSplittingCurrents))];
   Kirchhoff = Kirchhoff /. Join[RuleExitCurrentsIn, RuleEntryOut];
   {BM, KM} = 
    CoefficientArrays[Kirchhoff, 
     vars = Variables[Kirchhoff /. Equal -> Plus]];
   (*Print["The matrices B and K are: \n", MatrixForm /@ {-BM, KM}, 
    "\nThe order of the variables is \n", vars];*)
    cost = Function[currents,  MapThread[#1[#2]&, 
   	{KeyMap[Join[jvars, jtvars]][CostArgs] /@ vars, currents}]];
   {-BM, KM, cost, vars}
   ];

MFGPreprocessing[Eqs_] := 
  Module[{InitRules, RuleBalanceGatheringCurrents, EqEntryIn, 
    RuleEntryOut, RuleExitCurrentsIn, RuleExitValues, 
    EqValueAuxiliaryEdges, EqSwitchingByVertex, EqCompCon, 
    EqBalanceSplittingCurrents, EqCurrentCompCon, EqTransitionCompCon,
     EqPosJs, EqPosJts, ModuleVarsNames, ModulesVars, NewSystem}, 
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
   (*TODO Check switching costs for consistency*)(*First rules:
   entry currents*)
   InitRules = Association[Flatten[ToRules /@ EqEntryIn]];
   (*no exit at the entrances*)
   AssociateTo[InitRules, RuleEntryOut];
   (*no entrance at the exits*)
   AssociateTo[InitRules, RuleExitCurrentsIn];
   (*currents gathered from transition currents*)
   AssociateTo[InitRules, RuleBalanceGatheringCurrents];
   (*value function: exit costs*)
   AssociateTo[InitRules, RuleExitValues];
   EqSwitchingByVertex = 
    And @@ (Simplify /@ (Lookup[Eqs, "EqSwitchingByVertex", $Failed]
     /. InitRules));
   EqBalanceSplittingCurrents = 
    Lookup[Eqs, "EqBalanceSplittingCurrents", $Failed] /. InitRules;
   EqValueAuxiliaryEdges = 
    Lookup[Eqs, "EqValueAuxiliaryEdges", $Failed] /. InitRules;
    Rules = First @ Solve[EqBalanceSplittingCurrents
    && EqValueAuxiliaryEdges]//Quiet;
   InitRules = Join[InitRules/.Rules, Association@Rules];
   NewSystem = 
    MapThread[
     And, {{True, 
       EqPosJts && EqPosJs, 
       EqCurrentCompCon && EqTransitionCompCon && EqCompCon}, 
      Sys2Triple[EqSwitchingByVertex]}];
   {NewSystem, InitRules} = 
    TripleClean[{NewSystem, InitRules}];
   ModuleVarsNames = {"InitRules", "RuleBalanceGatheringCurrents", 
     "EqEntryIn", "RuleEntryOut", "RuleExitCurrentsIn", 
     "RuleExitValues", "EqValueAuxiliaryEdges", "EqSwitchingByVertex",
      "EqCompCon", "EqBalanceSplittingCurrents", "EqCurrentCompCon", 
     "EqTransitionCompCon", "EqPosJs", "EqPosJts", "NewSystem"};
   ModulesVars = {InitRules, RuleBalanceGatheringCurrents, EqEntryIn, 
      RuleEntryOut, RuleExitCurrentsIn, RuleExitValues, 
      EqValueAuxiliaryEdges, EqSwitchingByVertex, EqCompCon, 
      EqBalanceSplittingCurrents, EqCurrentCompCon, 
      EqTransitionCompCon, EqPosJs, EqPosJts, NewSystem} /. 
     InitRules;
       Join[Eqs, AssociationThread[ModuleVarsNames, ModulesVars]]];

CriticalCongestionSolver[Eqs_] := 
  Module[{PreEqs, Nlhs, EqCritical, InitRules, NewSystem,
  	vars, pickOne}, 
   PreEqs = MFGPreprocessing[Eqs];
   Nlhs = Lookup[PreEqs, "Nlhs", $Failed];
   InitRules = Lookup[PreEqs, "InitRules", $Failed];
   NewSystem = Lookup[PreEqs, "NewSystem", $Failed];
   (*Updated (with rules from preprocessing) edge equations*)
   EqCritical = (And @@ ((# == 0) & /@ Nlhs)) /. InitRules;
   Print[And[And@@NewSystem,EqCritical]];
   Print["FinalClean..."];
   {NewSystem, InitRules} = 
   FinalClean[{Join[{EqCritical}, Take[NewSystem, {2, 3}]], 
     InitRules}];
   (*Print["33: ", {NewSystem, InitRules}];*)
   (*NewSystem = NewReduce[And @@ NewSystem];
   Print["34: ", NewSystem];
   NewSystem = BooleanConvert[NewSystem, "CNF"];
   
   NewSystem = 
    Sys2Triple[BooleanConvert[Reduce[NewSystem, Reals], "CNF"]];
   Print["55: ",NewSystem];
   {NewSystem, InitRules} = TripleClean[{NewSystem, InitRules}];
   Print["66: ",NewSystem];
   {NewSystem, InitRules} = TripleClean[{Sys2Triple[Reduce[NewSystem,Reals]],InitRules}];*)
   NewSystem = And @@ NewSystem;
   Which[NewSystem === False, Print["There is no solution"],
   	NewSystem =!= True, 
   	NewSystem = Reduce[NewSystem, Reals];
   	{NewSystem, InitRules} = 
   FinalClean[{Sys2Triple[NewSystem], 
     InitRules}];
     NewSystem= And@@NewSystem;
   	Print["Multiple solutions: ", {NewSystem, InitRules}];
   	Print["\tPicking one..."];
	vars = Select[Join[Eqs["us"],Eqs["js"],Eqs["jts"]], Not[FreeQ[NewSystem, #]] &];
	(*Print[vars];*)
	(*Have to pick one so that all the currents have numerical values*)
	pickOne = Association@First@FindInstance[NewSystem && And @@ (# > 0 & /@vars), vars, Reals];
	InitRules = Expand/@Join[InitRules/.pickOne, pickOne];
   ];
   InitRules];

FinalStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, rules_}] := 
	{{EE, NN, OR}, rules}

FinalStep[{{EE_,NN_,OO_}, rules_}]:=
	Module[{NewSystem,newrules},
	{NewSystem,newrules} = TripleClean[{{EE,NN,OO}, rules}];
	Print["\tNewReduce...", {NewSystem}];
	NewSystem = NewReduce[And @@ NewSystem];
	Print["\tBooleanConvert...",NewSystem];
	NewSystem = BooleanConvert[NewSystem, "CNF"];
	NewSystem = Sys2Triple[NewSystem];
	{NewSystem, newrules}
	];

FinalClean[{{EE_, NN_, OR_}, rules_}] := 
 FixedPoint[FinalStep, {{EE, NN, OR}, rules}]

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

TripleStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, rules_}] := 
	{{EE, NN, OR}, rules}

TripleStep[{{EE_, NN_?TrueQ, OR_?TrueQ}, rules_Association}] := 
  Module[{newrules = {}}, newrules = First@Solve[EE] // Quiet;
   newrules = Join[rules /. newrules, Association@newrules];
   {{True, NN, OR}, Expand /@ newrules}];

TripleStep[{{EEs_, NNs_, ORs_}, rules_Association}] := 
  Module[{EE = EEs /. rules, 
  	NN = Simplify/@(NNs /. rules), 
  	OR = Simplify/@(ORs /. rules), 
    NNE, NNO, ORE, ORN, bool, newrules = {}, in, out}, 
    bool = EE && NN && OR;
   NN = Simplify[NN];
   {NNE, NN, NNO} = Sys2Triple[NN];
   {ORE, ORN, OR} = Sys2Triple[OR];
   EE = EE && NNE && ORE;
   If[EE =!= True && EE =!= False, 
    newrules = First@Solve[EE] // Quiet];
   newrules = Expand /@ Join[rules /. newrules, Association@newrules];
(*   in = Simplify/@((And[EEs,NNs,ORs]/.rules)&& (And @@ Equal @@@ Normal@rules));
   out = Simplify/@((And[EE,NN,OR]/.newrules)&& (And @@ Equal @@@ Normal@newrules));
  Print["Checking..."];
   Print[in];
  Print[out];
   *)
   
  (* in = Simplify/@(And[EEs,NNs,ORs]/.newrules);
   out = Simplify/@(And[EE,NN,OR]/.newrules);
   Print[in];
  Print[out];
  
   Print["Checking (replacing updated rules in both systems)..."];(**)   
    If[
   	Reduce[Implies[in, out],Reals]===True,
   	Print["Ok!"],
   	Print["Not Ok..."]
   ];	*)
   		
   		{{EE, NN, OR}, Expand /@ newrules}];


TripleClean[{{EE_, NN_, OR_}, rules_}] := 
 FixedPoint[TripleStep, {{EE, NN, OR}, rules}]


(*TODO review NewReduce and maybe make it take triple-systems*)

NewReduce[s_Or] := Reduce[s, Reals]

NewReduce[x_] := (Print["NewReduce[",x,"] = ", x];
  x)

(*TODO NewReduce:is there a problem here?*)

NewReduce[system_And] := 
 Module[{result, 
   groups = 
    GroupBy[List @@ system, Head[#] === Or || Head[#] === Equal &], 
   subgroups, sorted}, 
  subgroups = GroupBy[groups[True], Head[#] === Equal &];
  If[Head[subgroups[False]] === Missing, 
   result = 
    ZAnd[And @@ Lookup[groups, False, True], 
     And @@ Lookup[subgroups, True, True]], 
   sorted = 
    SortBy[And @@ Lookup[subgroups, False, True], 
     Simplify`SimplifyCount];
   result = 
    ZAnd[
     And @@ 
      Lookup[groups, False, 
       True], (And @@ Lookup[subgroups, True, True]) && sorted]];
  (*result=ZAnd[Select[system,!((Head[#]===Or)||(Head[#]===Equal))&],
  Select[system,((Head[#]===Or)||(Head[#]===Equal))&]];*)
  If[result =!= False, result = result // DeleteDuplicates];
  result]

ZAnd[_, False] := (False)

ZAnd[xp_, leq_LessEqual] := Simplify[xp && leq]

ZAnd[xp_, geq_GreaterEqual] := Simplify[xp && geq]

ZAnd[xp_, ineq_Inequality] := Simplify[xp && ineq]

ZAnd[False, _] := (False)

ZAnd[xp_, True] := xp

ZAnd[xp_, eq_Equal] := 
 With[{sol = Solve[eq]}, 
  If[sol === {}, 
   False, (xp /. First[sol]) && And @@ (First[sol] /. Rule -> Equal) //
     Simplify]]

ZAnd[xp_, 
  orxp_Or] := ((*Print[
   "Do we need removeduplicates?\n", (ZAnd[xp, #] & /@ orxp)];*)
  (ZAnd[xp, #] & /@ orxp) //RemoveDuplicates)

ZAnd[xp_, andxp_And] := 
 With[{fst = First[andxp], rst = Rest[andxp]}, 
  Which[
  	Head[fst] === Or, 
  	ReZAnd[xp, rst] /@ fst //RemoveDuplicates,
   	True, 
   	ReZAnd[xp, rst, fst]]]

(*Operator form of ReZAnd*)
ReZAnd[xp_, rst_] := ReZAnd[xp, rst, #] &

ReZAnd[xp_, rst_, fst_Equal] := 
 Module[{fsol = First@Solve@fst // Quiet}, 
  ZAnd[(xp /. fsol) && And @@ (fsol /. Rule -> Equal) // Simplify, 
   ReplaceSolution[rst, fsol]]]

ReZAnd[xp_, rst_, fst_And] := ReZAnd[xp, rst] /@ fst (*TODO never used???*)

ReZAnd[xp_, rst_, fst_] := 
 ZAnd[xp && fst, rst]

(*Inequality]:=ZAnd[xp&&fst,rst] ReZAnd[xp_,
  rst_,fst_LessEqual]:=ZAnd[xp&&fst,rst] ReZAnd[xp_,rst_,
  fst_GreaterEqual]:=ZAnd[xp&&fst,rst] ReZAnd[xp_,rst_,fst_Less]:=
  ZAnd[xp&&fst,rst] ReZAnd[xp_,rst_,fst_Greater*)
  
ReplaceSolution[rst_?BooleanQ, sol_] := rst

(*ReplaceSolution[False,sol_]:=False*)

ReplaceSolution[rst_, sol_] := Module[{newrst}, newrst = rst /. sol;
  If[Head[newrst] === And, Reduce[#, Reals] & /@ newrst, 
   Reduce[newrst, Reals]]]

(*ReplaceSolution[rst_,sol_]:=Simplify[rst/. sol]*)

RemoveDuplicates[xp_And]:=DeleteDuplicates[Sort[xp]];
RemoveDuplicates[xp_Or]:=DeleteDuplicates[Sort[xp]];
RemoveDuplicates[xp_]:=xp

(****************************************************************)
(****************************************************************)
(****************************************************************)
(****************************************************************)
(****************************************************************)
(****************************************************************)
(****************************************************************)
(****************************************************************)
(****************************************************************)



