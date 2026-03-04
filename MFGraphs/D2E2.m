(*Wolfram Language package*)

(* --- Switching cost consistency --- *)

ConsistentSwitchingCosts::usage =
"ConsistentSwitchingCosts[switchingcosts][{a,b,c}->S]
returns True if S, the cost of switching from the edge ab to cb, is smaller than any other combination,
such as, ab to bd and then from db to bc.
Returns the condition for this switching cost to satisfy the triangle inequality when S, and the other
switching costs too, does not have a numerical value.";

ConsistentSwitchingCosts[sc_][{a_, b_, c_} -> S_] :=
    Module[ {origin, bounds},
        origin = Cases[sc, HoldPattern[{a, b, _} -> _]];
        origin = DeleteCases[origin, {a, b, c} -> S];
        If[ origin =!= {},
            bounds = ((S <= Last[#] + Association[sc][{ Part[First[#], 3], b, c}]) & /@ origin);
            And @@ bounds,
            True
        ]
    ];

IsSwitchingCostConsistent::usage =
"IsSwitchingCostConsistent[List of switching costs] is True if all switching costs satisfy the triangle inequality. If some switching costs are symbolic, then it returns the consistency conditions."
IsSwitchingCostConsistent[switchingCosts_] :=
    And @@ Simplify[ConsistentSwitchingCosts[switchingCosts] /@ switchingCosts]

(* --- Graph helper functions --- *)

TransitionsAt[G_, k_] :=
    Insert[k,2] /@ Permutations[AdjacencyList[G, k], {2}]

AltFlowOp::usage =
"AltFlowOp[j][list] returns the alternative: j@@list ==0 || j@@Reverse@list ==0.";
AltFlowOp[j_][list_]:=
	j@@list==0||j@@Reverse@list==0;

FlowSplitting::usage =
"FlowSplitting[AT][UndirectedEdge[a, b]] returns the splitting that start with {a,b}.";
FlowSplitting[auxTriples_List][x_] := Select[auxTriples, MatchQ[#,{Sequence@@x,__}]&] ;

FlowGathering::usage=
"FlowGathering[auxTriples_List][x_] returns the triples that end with x.";
FlowGathering[auxTriples_List][x_] := Select[auxTriples, MatchQ[#,{__,Sequence@@x}]&]

IncomingEdges[FG_][k_] :=
    {k, #1} & /@ IncidenceList[FG, k];

IneqSwitch::usage =
"IneqSwitch[u, switchingCosts][{v,e1,e2}] returns the optimality condition at the vertex v related to switching from e1 to e2. Namely,
u[{v, e1}] <= u[{v, e2}] + switchingCosts[{v, e1, e2}]"
IneqSwitch[u_,Switching_Association][r_,i_,w_] := u[r,i] <= u[w,i]+Switching[{r,i,w}];

AltSwitch::usage =
"AltSwitch[jt,u,switchingCosts][{v,e1,e2}] returns the complementarity condition:
(jt[{v, e1, e2}] == 0) || (u[{v, e2}] - u[{v, e1}] + switchingCosts[{v, e1, e2}] == 0)"
AltSwitch[j_, u_, Switching_][r_, i_, w_] :=
    (j[r,i,w] == 0) ||
    (u[r,i] == u[w,i] + Switching[{r,i,w}]);

(* --- Data2Equations: main converter --- *)

Data2Equations::usage = "Data2Equations[Data] returns the equations, \
inequalities, and alternatives associated to the Data. "
Data2Equations[Data_Association] :=
    Module[ {verticesList, adjacencyMatrix, entryVerticesFlows, exitVerticesCosts,
    	switchingCosts, graph, entryVertices, auxEntryVertices,
      exitVertices, auxExitVertices, entryEdges, exitEdges,
      auxiliaryGraph, auxEdgeList, edgeList, auxVerticesList, auxPairs, auxTriples,
      EntryDataAssociation, ExitCosts, js, us, jts,
      SignedFlows, SwitchingCosts, IneqJs, IneqJts,
      AltFlows, AltTransitionFlows, splittingPairs,
      EqBalanceSplittingFlows, BalanceSplittingFlows,
      NoDeadStarts, RuleBalanceGatheringFlows,
      BalanceGatheringFlows, EqBalanceGatheringFlows, EqEntryIn, RuleEntryValues,
      RuleEntryOut, RuleExitFlowsIn, RuleExitValues,
      EqValueAuxiliaryEdges, IneqSwitchingByVertex,
      AltOptCond, Nlhs, ModuleVars, ModuleVarsNames,
      Nrhs, costpluscurrents, EqGeneral,
      inAuxEntryPairs, outAuxEntryPairs, inAuxExitPairs,
      outAuxExitPairs, pairs, halfPairs,
      consistentCosts},
        verticesList = Lookup[Data, "Vertices List", {}];
        adjacencyMatrix = Lookup[Data, "Adjacency Matrix", {}];
        entryVerticesFlows = Lookup[Data, "Entrance Vertices and Flows", {}];
        exitVerticesCosts = Lookup[Data, "Exit Vertices and Terminal Costs", {}];
        switchingCosts = Lookup[Data, "Switching Costs", {}];

        (* Graph construction *)
        graph = AdjacencyGraph[verticesList, adjacencyMatrix, VertexLabels -> "Name", DirectedEdges -> False];
        entryVertices = First /@ entryVerticesFlows;
        exitVertices = First /@ exitVerticesCosts;

        (* Define auxiliary vertices for the entry and exit vertices *)
        auxEntryVertices = Symbol["en" <> ToString[#]] &/@ entryVertices;
        auxExitVertices = Symbol["ex" <> ToString[#]] &/@ exitVertices;

        entryEdges = MapThread[UndirectedEdge, {auxEntryVertices, entryVertices}];
        exitEdges = MapThread[UndirectedEdge, {exitVertices, auxExitVertices}];
        auxiliaryGraph = EdgeAdd[graph, Join[entryEdges, exitEdges]];
        auxEdgeList = EdgeList[auxiliaryGraph];
        edgeList = EdgeList[graph];
        auxVerticesList = VertexList[auxiliaryGraph];

        (* Arguments for the relevant functions: flows, values, and transition flows *)
        halfPairs = List@@@edgeList;
        inAuxEntryPairs = List@@@entryEdges;
        outAuxExitPairs = List@@@exitEdges;

        inAuxExitPairs = Reverse/@outAuxExitPairs;
        outAuxEntryPairs = Reverse/@inAuxEntryPairs;
        pairs = Join[halfPairs, Reverse/@halfPairs];
        auxPairs = Join[inAuxEntryPairs, outAuxEntryPairs, inAuxExitPairs, outAuxExitPairs, pairs];
        (* Insert the vertex in pairs of adjacent vertices *)
        auxTriples = Flatten[Insert[#, 2] /@ Permutations[AdjacencyList[auxiliaryGraph, #], {2}]& /@ auxVerticesList, 1];

        (* Prepare boundary data *)
        EntryDataAssociation = RoundValues@AssociationThread[inAuxEntryPairs, Last /@ entryVerticesFlows];
        ExitCosts = AssociationThread[auxExitVertices, Last /@ exitVerticesCosts];

        (* Variables *)
        js = j[Sequence @@ #] & /@ auxPairs;
        jts = j[Sequence @@ #] & /@ auxTriples;
        us = u[Sequence @@ #] & /@ auxPairs;

        (* Signed flows for "half" of the variables *)
        SignedFlows = AssociationMap[j@@# - j@@Reverse@# &, Join[inAuxEntryPairs, outAuxExitPairs, halfPairs]];
        SwitchingCosts = AssociationMap[0 &, auxTriples];
        If[switchingCosts =!= {},
        	AssociateTo[SwitchingCosts, AssociationThread[Most /@ switchingCosts, Last /@ switchingCosts]]
        ];
        consistentCosts = IsSwitchingCostConsistent[Normal@SwitchingCosts];
        Which[consistentCosts === False,
            Print["Switching costs are inconsistent!"];
            Return[ $Failed, Module],
         consistentCosts =!= True,
         MFGPrint["Switching costs conditions are ", consistentCosts]
        ];

        IneqJs = And @@ (# >= 0 & /@ Join[js]);
        IneqJts = And @@ (# >= 0 & /@ Join[jts]);
        AltFlows = And@@(AltFlowOp[j]/@auxPairs);
        AltTransitionFlows = And@@(AltFlowOp[j]/@auxTriples);
        splittingPairs = Join[inAuxEntryPairs, inAuxExitPairs, pairs];
        BalanceSplittingFlows = (j@@# - Total[j@@@FlowSplitting[auxTriples][#]])& /@ splittingPairs;
        EqBalanceSplittingFlows = Simplify /@ (And @@ ((# == 0) & /@ BalanceSplittingFlows));
        NoDeadStarts = Join[outAuxEntryPairs, outAuxExitPairs, pairs];
        RuleBalanceGatheringFlows = Association[(j@@# -> Total[j@@@ FlowGathering[auxTriples][#]]) & /@ NoDeadStarts];

        (* Equations for the exit currents at the entry vertices *)
        BalanceGatheringFlows = ((-j@@# + Total[j@@@ FlowGathering[auxTriples][#]]) & /@ NoDeadStarts);
        EqBalanceGatheringFlows = Simplify /@ (And @@ (# == 0 & /@ BalanceGatheringFlows));

        (* Incoming currents *)
        EqEntryIn = (j@@# == EntryDataAssociation[#]) & /@  inAuxEntryPairs;

        (* Outgoing flows at entrances and incoming flows at exits are zero *)
        RuleEntryOut = Association[(j@@# -> 0) & /@ outAuxEntryPairs];
        RuleExitFlowsIn = Association[(j@@# -> 0) & /@ inAuxExitPairs];

		(* Exit values at exit vertices *)
        RuleExitValues = AssociationThread[u@@@(Reverse/@outAuxExitPairs),Last /@ exitVerticesCosts];
        RuleExitValues = Join[RuleExitValues, AssociationThread[u@@@outAuxExitPairs, Last /@ exitVerticesCosts]];
        RuleEntryValues = AssociationThread[u@@@outAuxEntryPairs, u@@@inAuxEntryPairs];
        EqValueAuxiliaryEdges = And @@ ((u@@# == u@@Reverse[#]) & /@ Join[inAuxEntryPairs, outAuxExitPairs]);
        IneqSwitchingByVertex = IneqSwitch[u, SwitchingCosts] @@@ TransitionsAt[auxiliaryGraph, #] & /@ verticesList;
        IneqSwitchingByVertex = And @@@ IneqSwitchingByVertex;
        AltOptCond = And @@ AltSwitch[j, u, SwitchingCosts] @@@ auxTriples;
        Nlhs = Flatten[u@@# - u@@Reverse@# + SignedFlows[#] & /@ halfPairs];
        Nrhs = Flatten[SignedFlows[#] - Sign[SignedFlows[#]] Cost[SignedFlows[#], #] & /@ halfPairs];

       	(* Cost-plus-currents for the general (non-critical) case *)
        costpluscurrents = Table[Symbol["cpc" <> ToString[k]], {k, 1, Length@edgeList}];
        EqGeneral = And @@ (MapThread[Equal, {Nlhs, costpluscurrents}]);
        costpluscurrents = AssociationThread[costpluscurrents, Nrhs];

        (* Build output association *)
        ModuleVars = {graph, pairs,
          entryVertices, auxEntryVertices, exitVertices, auxExitVertices,
          entryEdges, exitEdges, auxiliaryGraph, auxEdgeList, edgeList, auxVerticesList,
          auxTriples, EntryDataAssociation,
          ExitCosts, js, us, jts, SignedFlows,
          SwitchingCosts, IneqJs, IneqJts, AltFlows,
          AltTransitionFlows, splittingPairs, EqBalanceSplittingFlows,
          BalanceSplittingFlows, NoDeadStarts,
          RuleBalanceGatheringFlows, BalanceGatheringFlows,
          EqBalanceGatheringFlows, EqEntryIn, RuleEntryOut, RuleEntryValues,
          RuleExitFlowsIn, RuleExitValues, EqValueAuxiliaryEdges,
          IneqSwitchingByVertex, AltOptCond, Nlhs,
          Nrhs, costpluscurrents, EqGeneral};
        ModuleVarsNames = {"graph", "pairs",
          "entryVertices", "auxEntryVertices", "exitVertices",
          "auxExitVertices", "entryEdges", "exitEdges", "auxiliaryGraph",
           "auxEdgeList", "edgeList", "auxVerticesList", "auxTriples",
           "EntryDataAssociation", "ExitCosts", "js",
          "us", "jts", "SignedFlows",
          "SwitchingCosts", "IneqJs", "IneqJts", "AltFlows",
          "AltTransitionFlows", "splittingPairs",
          "EqBalanceSplittingFlows", "BalanceSplittingFlows",
          "NoDeadStarts", "RuleBalanceGatheringFlows",
          "BalanceGatheringFlows", "EqBalanceGatheringFlows",
          "EqEntryIn", "RuleEntryOut", "RuleEntryValues", "RuleExitFlowsIn",
          "RuleExitValues", "EqValueAuxiliaryEdges",
           "IneqSwitchingByVertex", "AltOptCond", "Nlhs",
          "Nrhs", "costpluscurrents", "EqGeneral"};
        Join[Data, AssociationThread[ModuleVarsNames, ModuleVars]]
    ];

(* --- Utility --- *)

NumberVectorQ::usage =
"NumberVectorQ[j] returns True if the vector j is numeric."
NumberVectorQ[j_] :=
    And @@ (NumberQ /@ j);

RoundValues[x_?NumberQ] :=
    Round[x, 10^-10]

RoundValues[Rule[a_, b_]] :=
    Rule[a, RoundValues[b]]

RoundValues[x_List] :=
    RoundValues /@ x

RoundValues[x_Association] :=
    RoundValues /@ x

(* --- GetKirchhoffMatrix --- *)

GetKirchhoffMatrix::usage = "GetKirchhoffMatrix[d2e] returns the \
entry current vector, Kirchhoff matrix, (critical congestion) cost \
function, and the variables in the order corresponding to the Kirchhoff matrix."
GetKirchhoffMatrix[Eqs_] :=
    Module[ {Kirchhoff, EqEntryIn = Lookup[Eqs, "EqEntryIn", True],
      BalanceGatheringFlows =
       Lookup[Eqs, "BalanceGatheringFlows", {}],
      BalanceSplittingFlows =
       Lookup[Eqs, "BalanceSplittingFlows", {}],
      RuleExitFlowsIn = Lookup[Eqs, "RuleExitFlowsIn", {}],
      RuleEntryOut = Lookup[Eqs, "RuleEntryOut", {}], BM, KM, vars,
      CostArgs = Lookup[Eqs, "CostArgs", <||>],
      jvars = Lookup[Eqs, "jvars", {}],
      jtvars = Lookup[Eqs, "jtvars", {}], cost, CCost},
        Kirchhoff = Join[EqEntryIn, (# == 0 & /@ (BalanceGatheringFlows + BalanceSplittingFlows))];
        Kirchhoff = Kirchhoff /. Join[RuleExitFlowsIn, RuleEntryOut];
        vars = Select[Values@Join[jvars,jtvars],MemberQ[Variables[Kirchhoff /. Equal -> List],#]&];
        {BM, KM} = CoefficientArrays[Kirchhoff, vars];
        cost = AssociationThread[vars, KeyMap[Join[jvars,jtvars]][CostArgs]/@vars];
        CCost = cost /@ vars /. MapThread[Rule, {vars, #}] &;
        {-BM, KM, CCost, vars}
    ];

(* --- MFGPreprocessing --- *)

MFGPreprocessing::usage =
"MFGPreprocessing[Eqs] returns the association Eqs with the preliminary solution \"InitRules\" and corresponding 'reduced' \"NewSystem\"."
MFGPreprocessing[Eqs_] :=
    Module[ {InitRules, RuleBalanceGatheringFlows, EqEntryIn,
      RuleEntryOut, RuleExitFlowsIn, RuleExitValues,
      EqValueAuxiliaryEdges, IneqSwitchingByVertex, AltOptCond,
      EqBalanceSplittingFlows, AltFlows, AltTransitionFlows,
      IneqJs, IneqJts, ModuleVarsNames, ModulesVars, NewSystem,
      Rules, EqGeneral, RuleEntryValues, temp, time},
        RuleBalanceGatheringFlows = Lookup[Eqs, "RuleBalanceGatheringFlows", $Failed];
        EqEntryIn = Lookup[Eqs, "EqEntryIn", $Failed];
        RuleEntryOut = Lookup[Eqs, "RuleEntryOut", $Failed];
        RuleExitFlowsIn = Lookup[Eqs, "RuleExitFlowsIn", $Failed];
        RuleExitValues = Lookup[Eqs, "RuleExitValues", $Failed];
        RuleEntryValues = Lookup[Eqs, "RuleEntryValues", $Failed];
        AltOptCond = Lookup[Eqs, "AltOptCond", $Failed];
        AltFlows = Lookup[Eqs, "AltFlows", $Failed];
        AltTransitionFlows = Lookup[Eqs, "AltTransitionFlows", $Failed];
        IneqJs = Lookup[Eqs, "IneqJs", $Failed];
        IneqJts = Lookup[Eqs, "IneqJts", $Failed];
        EqGeneral = Lookup[Eqs, "EqGeneral", $Failed];
        IneqSwitchingByVertex = Lookup[Eqs, "IneqSwitchingByVertex", $Failed];
        EqBalanceSplittingFlows = Lookup[Eqs, "EqBalanceSplittingFlows", $Failed];
        EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", $Failed];

        MFGPrint[Eqs["auxiliaryGraph"]];
        (* First rules: entry currents *)
        InitRules = Association[Flatten[ToRules /@ EqEntryIn]];
        AssociateTo[InitRules, RuleEntryOut];
        AssociateTo[InitRules, RuleExitFlowsIn];
        AssociateTo[InitRules, RuleEntryValues];
        AssociateTo[InitRules, RuleBalanceGatheringFlows];
        AssociateTo[InitRules, RuleExitValues];

        temp = MFGPrintTemporary["Simplifying inequalities involving Switching Costs...", IneqSwitchingByVertex];
        {time, IneqSwitchingByVertex} = AbsoluteTiming[And @@ (Simplify /@ ( IneqSwitchingByVertex /. InitRules))];
        NotebookDelete[temp];
        MFGPrint["Switching costs simplified in ", time, " seconds."];
        EqBalanceSplittingFlows = EqBalanceSplittingFlows /. InitRules;
        EqValueAuxiliaryEdges = EqValueAuxiliaryEdges /. InitRules;
        temp = MFGPrintTemporary["Solving some balance equations: ", EqBalanceSplittingFlows && EqValueAuxiliaryEdges];
        {time, Rules} = AbsoluteTiming[First@Solve[EqBalanceSplittingFlows && EqValueAuxiliaryEdges] // Quiet];
        InitRules = Join[InitRules /. Rules, Association@Rules];
        NotebookDelete[temp];
        MFGPrint["Balance equations solved in ", time, " seconds."];
        NewSystem =
         MapThread[
          And, {{True, IneqJts && IneqJs, AltFlows && AltTransitionFlows && AltOptCond},
           Sys2Triple[IneqSwitchingByVertex]}
         ];
        temp = MFGPrintTemporary["TripleClean will work on\n", NewSystem, "\nwith: ", InitRules];
        {time,{NewSystem, InitRules}} = AbsoluteTiming@TripleClean[{NewSystem, InitRules}];
        NotebookDelete[temp];
        MFGPrint["TripleClean took ", time, " seconds."];
        NewSystem[[1]] = EqGeneral;
        temp = MFGPrintTemporary["TripleClean again (with some more equalities)..."];
        {time,{NewSystem, InitRules}} = AbsoluteTiming@TripleClean[{NewSystem, InitRules}];
        NotebookDelete[temp];
        MFGPrint["TripleClean again (with some more equalities) took ", time, " seconds."];

        ModuleVarsNames = {"InitRules", "NewSystem"};
        ModulesVars = {InitRules, NewSystem};
        Join[Eqs, AssociationThread[ModuleVarsNames, ModulesVars]]
    ];

(* --- CriticalCongestionSolver --- *)

CriticalCongestionSolver::usage =
  "CriticalCongestionSolver[Eqs] returns Eqs with an association \"AssoCritical\" with rules to
the solution to the critical congestion case";

CriticalCongestionSolver[$Failed] := $Failed

CriticalCongestionSolver[Eqs_] :=
    Module[ {PreEqs, js, AssoCritical, time, temp},
    	ClearSolveCache[];
    	If[KeyExistsQ[Eqs,"InitRules"],
    		PreEqs = Eqs,
    		temp = MFGPrintTemporary["Preprocessing..."];
        	{time, PreEqs} = AbsoluteTiming@MFGPreprocessing[Eqs];
        	NotebookDelete[temp];
        	MFGPrint["Preprocessing took ", time, " seconds to terminate."];
        ];
        js = Lookup[PreEqs, "js", $Failed];
        AssoCritical = MFGSystemSolver[PreEqs][AssociationThread[js, 0 js]];
        Join[PreEqs, Association["AssoCritical"-> AssoCritical]]
    ];

(* --- MFGSystemSolver --- *)

MFGSystemSolver::usage =
"MFGSystemSolver[Eqs][edgeEquations] returns the
association with rules to the solution";
MFGSystemSolver[Eqs_][approxJs_] :=
    Module[ {NewSystem, InitRules, pickOne, vars, System, Ncpc,
        costpluscurrents, us, js, jts, jjtsR, usR, time, temp,
        ineqsByTransition},
        ClearSolveCache[];
        us = Lookup[Eqs, "us", $Failed];
        js = Lookup[Eqs, "js", $Failed];
        jts = Lookup[Eqs, "jts", $Failed];
        InitRules = Lookup[Eqs, "InitRules", $Failed];
        NewSystem = Lookup[Eqs, "NewSystem", $Failed];
        costpluscurrents = Lookup[Eqs, "costpluscurrents", $Failed];
        {time, Ncpc} = AbsoluteTiming[RoundValues @ (Expand/@(costpluscurrents /. approxJs))];
        MFGPrint["MFGSS: Calculated the cost plus currents for the flow in ", time, " seconds."];
        InitRules = Expand/@(InitRules /. Ncpc);
        NewSystem = NewSystem /. Ncpc;

        If[And@@((#===False)&/@NewSystem),
        	Print["MFGSS: There is no solution!"];
            Return[Null,Module]
        ];

        (* Retrieve some equalities from the inequalities: group by transition flow *)
        temp = MFGPrintTemporary["MFGSS: Selecting inequalities by transition flow..."];
        {time, ineqsByTransition} = AbsoluteTiming[Select[NewSystem[[2]], Function[exp, !FreeQ[#][exp]]]&/@jts];
        NotebookDelete[temp];
        MFGPrint["MFGSS: Selecting inequalities by transition flow took ", time, " seconds. ", Length[ineqsByTransition]];
        temp = MFGPrintTemporary["MFGSS: Simplifying inequalities by transition flow..."];
        {time, ineqsByTransition}= AbsoluteTiming[RemoveDuplicates[Simplify/@ineqsByTransition]];
        NotebookDelete[temp];
        MFGPrint["MFGSS: Simplifying inequalities by transition flow took ", time, " seconds. "];

		NewSystem[[2]]=(And@@ineqsByTransition);
		NewSystem = Sys2Triple[And@@NewSystem];
		{NewSystem,InitRules} = TripleClean[{NewSystem,InitRules}];

		{NewSystem, InitRules} = FinalStep[{NewSystem, InitRules}];

        System = BooleanConvert[NewSystem[[3]]]&&NewSystem[[2]];
        MFGPrint["BooleanConvert done. Simplifying..."];
        System = Simplify[And@@NewSystem];
        MFGPrint["Simplifying done. TripleClean..."];
        {System, InitRules} = TripleClean[{Sys2Triple[System], InitRules}];
        System = And@@System;
        Which[System === False,
            Print["MFGSS: There is no solution!"];
            Return[Null,Module],
            System =!= True,
            MFGPrint["MFGSS: (Possibly) Multiple solutions:\n", System];
            jjtsR = Join[js, jts];
            jjtsR = Intersection[Flatten[List@@@((Join[us, js] /. InitRules)/.Times->Plus)], jjtsR];
        	usR = Intersection[Flatten[List@@@((Join[us, js] /. InitRules)/.Times->Plus)], us];
            vars = Join[usR, jjtsR];

            (* Pick one solution so that all the currents have numerical values *)
            If[Length[vars] > 0,
              pickOne = FindInstance[System && And @@ ((# > 0 )& /@ jjtsR), vars, Reals];
              If[pickOne === {},
              	pickOne = FindInstance[System && And @@ ((# >= 0 )& /@ jjtsR), vars, Reals]
              ];
              If[pickOne =!= {},
                pickOne = Association @ First @ pickOne;
                MFGPrint["MFGSS: Picked one solution: ", pickOne];
                InitRules = Expand /@ Join[InitRules /. pickOne, pickOne],
                MFGPrint["MFGSS: No feasible solution found"];
                Return[Null, Module]
              ],
              (* All variables already solved — nothing to pick *)
              MFGPrint["MFGSS: All variables already determined in InitRules"]
            ],
            True,
         	MFGPrint["MFGSS: System is ", System]
         ];
        InitRules = Join[KeyTake[InitRules, us], KeyTake[InitRules, js], KeyTake[InitRules, jts]];
        InitRules
    ];

(* --- FinalStep --- *)

FinalStep::usage =
"FinalStep[{EE,NN,OR}, rules] takes a grouped system and some Association of rules (a partial solution). It returns the result of applying DNFReduce.";

FinalStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, rules_}] :=
    {{EE, NN, OR}, rules}

FinalStep[{{EE_, NN_, OO_}, rules_}] :=
    Module[ {NewSystem, newrules, sorted = True, time, temp},
        {NewSystem, newrules} = TripleClean[{{EE, NN, OO}, rules}];
        If[NewSystem[[3]] === True,
            Return[{NewSystem, newrules}, Module],
            MFGPrint["Final: ", Length/@NewSystem];
            sorted = RemoveDuplicates[NewSystem[[3]]];
        ];
        temp = MFGPrintTemporary["Final: Iterative DNF conversion on " , Length[sorted]," disjunctions..."];
        {time, NewSystem} = AbsoluteTiming[DNFReduce[And @@ Most[NewSystem], sorted]];
        NotebookDelete[temp];
        MFGPrint["Final: Iterative DNF conversion on " , Length[sorted]," disjunctions took ", time, " seconds to terminate."];
        If[Head[NewSystem] === Or,
	        MFGPrint["Final: Reducing (each alternative)... ", Length[NewSystem] ," of them."];
        	{time, NewSystem} = AbsoluteTiming[Reduce[#,Reals]&/@NewSystem];
        	NewSystem = Simplify[NewSystem],
        	{time,NewSystem} = AbsoluteTiming@Reduce[NewSystem, Reals]
        ];
        MFGPrint["Final: Reducing (each alternative), returned ", Length@NewSystem," of them, ","took ", time, " seconds to terminate."];
        NewSystem = BooleanConvert@NewSystem;
        NewSystem = Sys2Triple[NewSystem];
        {NewSystem, newrules} = TripleClean[{NewSystem, newrules}];
        MFGPrint["Now: ", TimeObject[Now], " The new rules are: ", newrules, Length/@NewSystem];
        {NewSystem, newrules}
    ];

(* --- Sys2Triple: decompose into equalities, inequalities, alternatives --- *)

Sys2Triple::usage =
"Sys2Triple[sys] returns a triple {equalities, inequalities, alternatives} from sys.
The input should be a system of equations, inequalities and (simple) alternatives."

Sys2Triple[True] = Table[True, 3]

Sys2Triple[False] = Table[False, 3]

Sys2Triple[system_] :=
    Which[
     Head[system] === And,
      Module[ {groups, EE, OR, NN},
          groups = GroupBy[List @@ system, Head[#] === Equal &];
          EE = And @@ Lookup[groups, True, {}];
          groups = GroupBy[Lookup[groups, False, {}], Head[#] === Or &];
          OR = And @@ Lookup[groups, True, True];
          NN = And @@ Lookup[groups, False, True];
          {EE, NN, OR}
      ],
     Head[system] === Equal,
      {system, True, True},
     Head[system] === Or,
         {True, True, system},
     True,
         {True, system, True}];

(* --- TripleStep / TripleClean: fixed-point simplification --- *)

TripleStep::usage =
"TripleStep[{{EE,NN,OR},Rules}] returns {{NewEE, NewNN, NewOR}, NewRules}, where NewRules contain the solutions to all the equalities found in the system after replacing Rules in {EE,NN,OR}."

TripleStep[{{EEs_, NNs_, ORs_}, rules_List}] :=
    TripleStep[{{EEs, NNs, ORs}, Association@rules}]

TripleStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, rules_Association}] :=
    {{EE, NN, OR}, rules}

TripleStep[{{EE_, NN_?TrueQ, OR_?TrueQ}, rules_Association}] :=
    Module[ {newrules = {}},
        newrules = First@Solve[EE /. rules] // Quiet;
        newrules = Join[rules /. newrules, Association@newrules];
        {{True, NN, OR}, Expand /@ newrules}
    ];

TripleStep[{{True, NNs_, ORs_}, rules_Association}] :=
	{{True, NNs, ORs}, rules}

TripleStep[{{EEs_, NNs_, ORs_}, rules_Association}] :=
	Module[{EE = EEs /. rules, NN, OR, NNE, NNO, ORE, ORN, newrules ={}},
	If[ EE =!= True && EE =!= False,
    	newrules = First@Solve[EE] // Quiet
    ];
    newrules = Expand /@ Join[rules /. newrules, Association@newrules];
    NN = Expand /@ (NNs /. newrules);
    OR = Expand /@ (ORs /. newrules);
    {NNE, NN, NNO} = Sys2Triple[NN];
    {ORE, ORN, OR} = Sys2Triple[OR];
    EE = NNE && ORE;
    {{EE, NN, OR}, newrules}
    ];

TripleClean::usage =
"TripleClean[{{EE,NN,OR},Rules}] composes TripleStep until it reaches a fixed point, that is, {{True,NewNN,NewOR},NewRules} such that replacement of NewRules in NewNN and NewOR do not produce equalities."
TripleClean[{{EE_, NN_, OR_}, rules_}] := FixedPoint[TripleStep, {{EE, NN, OR}, rules}];
