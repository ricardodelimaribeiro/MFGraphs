(*Wolfram Language package*)
ConsistentSwithingCosts::usage = 
"ConsistentSwithingCosts[switchingcosts][{a,b,c}->S]  
returns True if S, the cost of switching from the edge ab to cb, is smaller than any other combination, such as, ab to bd and then from db to bc.
returns the condition for this switching cost to satisfy the triangle inequality when S, and the other switching costs too, does not have a numerical value.";

ConsistentSwithingCosts[sc_][{a_, b_, c_} -> S_] :=
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
    And @@ Simplify[ConsistentSwithingCosts[switchingCosts] /@ switchingCosts]

TransitionsAt[G_, k_] :=
    Insert[k,2] /@ Permutations[AdjacencyList[G, k], {2}]


AltFlowOp::usage = 
"AltFlowOp[j][list] returns the alternative: j@@list ==0 || j@@Reverse@list ==0.";     
AltFlowOp[j_][list_]:=
	j@@list==0||j@@Reverse@list==0;
	
(*Clear[FlowSplitting];*)
(*TODO: resolve kirchhoff equations: *)
FlowSplitting::usage = 
"FlowSplitting[AT][UndirectedEdge[a, b]] returns the splitting that start with {a,b}.";
FlowSplitting[auxTriples_List][x_] := Select[auxTriples, MatchQ[#,{Sequence@@x,__}]&] ;


FlowGathering::usage=
"FlowGathering";
FlowGathering[auxTriples_List][x_] := Select[auxTriples, MatchQ[#,{__,Sequence@@x}]&]

IncomingEdges[FG_][k_] :=
    {k, #1} & /@ IncidenceList[FG, k];

OutgoingEdges[FG_][k_] :=
    OtherWay /@ ({k, #} & /@ IncidenceList[FG, k]);

IneqSwitch::usage = 
"IneqSwitch[u, switchingCosts][{v,e1,e2}] returns the optimality condition at the vertex v related to switching form e1 to e2. Namely, 
u[{v, e1}] <= u[{v, e2}] + switchingCosts[{v, e1, e2}]"
IneqSwitch[u_,Switching_Association][r_,i_,w_] := u[r,i] <= u[w,i]+Switching[{r,i,w}];

AltSwitch::usage =
"AltSwitch[jt,u,switchingCosts][{v,e1,e2}] returns the complementarity condition:
(jt[{v, e1, e2}] == 0) || (u[{v, e2}] - u[{v, e1}] + switchingCosts[{v, e1, e2}] == 0)"
AltSwitch[j_, u_, Switching_][r_, i_, w_] :=
    (j[r,i,w] == 0) || 
    (u[r,i] == u[w,i] + Switching[{r,i,w}]);

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
      EqValueAuxiliaryEdges, OutRules, InRules, IneqSwitchingByVertex, 
      AltOptCond, Nlhs, ModuleVars, ModuleVarsNames, 
      CostArgs, Nrhs, 
      consistentCosts, costpluscurrents, EqGeneral, 
      inAuxEntryPairs, outAuxEntryPairs, inAuxExitPairs, 
      outAuxExitPairs, pairs, halfPairs},
        verticesList = Lookup[Data, "Vertices List", {}];
        adjacencyMatrix = Lookup[Data, "Adjacency Matrix", {}];
        entryVerticesFlows = Lookup[Data, "Entrance Vertices and Flows", {}];
        exitVerticesCosts = Lookup[Data, "Exit Vertices and Terminal Costs", {}];
        switchingCosts = Lookup[Data, "Switching Costs", {}];
        (*Print[entryVerticesFlows];*)
        (****Graph stuff****)
        graph = AdjacencyGraph[verticesList, adjacencyMatrix, VertexLabels -> "Name", DirectedEdges -> False];
        entryVertices = First /@ entryVerticesFlows;
        exitVertices = First /@ exitVerticesCosts;
        
        (*Define auxiliary vertices for the entry and exit vertices*)
        auxEntryVertices = Symbol["en" <> ToString[#]] &/@ entryVertices;
        auxExitVertices = Symbol["ex" <> ToString[#]] &/@ exitVertices;
        
        entryEdges = MapThread[UndirectedEdge, {auxEntryVertices, entryVertices}];
        exitEdges = MapThread[UndirectedEdge, {exitVertices, auxExitVertices}];
        auxiliaryGraph = EdgeAdd[graph, Join[entryEdges, exitEdges]];
        (*Print[auxiliaryGraph];*)
        auxEdgeList = EdgeList[auxiliaryGraph];
        edgeList = EdgeList[graph];
        auxVerticesList = VertexList[auxiliaryGraph];
        
        (*arguments for the relevant functions on the graph: flows, values, and transition flows.*)
        halfPairs = List@@@edgeList;
        inAuxEntryPairs = List@@@entryEdges;
        outAuxExitPairs = List@@@exitEdges;
        
        inAuxExitPairs = Reverse/@outAuxExitPairs;
        outAuxEntryPairs = Reverse/@inAuxEntryPairs;
        pairs = Join[halfPairs, Reverse/@halfPairs];
        auxPairs = Join[inAuxEntryPairs, outAuxEntryPairs, inAuxExitPairs, outAuxExitPairs, pairs];
        (*Insert the vertex in pairs of adjacent vertices*)
        auxTriples = Flatten[Insert[#, 2] /@ Permutations[AdjacencyList[auxiliaryGraph, #], {2}]& /@ auxVerticesList, 1];
        
        (*Prepare boundary data*)
        EntryDataAssociation = RoundValues@AssociationThread[inAuxEntryPairs, Last /@ entryVerticesFlows];
        ExitCosts = AssociationThread[auxExitVertices, Last /@ exitVerticesCosts];
        (*Print[auxExitVertices,ExitCosts];*)
        
        (*variables*)
        js = j[Sequence @@ #] & /@ auxPairs;
        jts = j[Sequence @@ #] & /@ auxTriples;
        us = u[Sequence @@ #] & /@ auxPairs;
        
        (*These are the signed flows for "half" of the variables.*)
        SignedFlows = AssociationMap[j@@# - j@@Reverse@# &, Join[inAuxEntryPairs, outAuxExitPairs, halfPairs]];
        (*Print["SF: ",SignedFlows];*)
        SwitchingCosts = AssociationMap[0 &, auxTriples];
        If[switchingCosts =!= {}, 
        	AssociateTo[SwitchingCosts, AssociationThread[Most /@ switchingCosts, Last /@ switchingCosts]]
        ];
        (*Print["switchingCosts: ", SwitchingCosts];*)
        consistentCosts = IsSwitchingCostConsistent[Normal@SwitchingCosts];
        Which[consistentCosts === False, 
            Print["Switching costs are inconsistent!"];
            Return[ $Failed, Module], 
         consistentCosts =!= True, 
         Print["Switching costs conditions are ", consistentCosts]
        ];
        
        
        IneqJs = And @@ (# >= 0 & /@ Join[js]);(*Inequality*)
        IneqJts = And @@ (# >= 0 & /@ Join[jts]);(*Inequality*)
        AltFlows = And@@(AltFlowOp[j]/@auxPairs);
        AltTransitionFlows = And@@(AltFlowOp[j]/@auxTriples);
        (*Print[IneqJs,IneqJts,AltFlows,AltTransitionFlows];*)
        splittingPairs = Join[inAuxEntryPairs, inAuxExitPairs, pairs];
        BalanceSplittingFlows = (j@@# - Total[j@@@FlowSplitting[auxTriples][#]])& /@ splittingPairs;
        EqBalanceSplittingFlows = Simplify /@ (And @@ ((# == 0) & /@ BalanceSplittingFlows));
        (*Print["nde: ",splittingPairs, "\n",BalanceSplittingFlows,"\nEBSF: ",EqBalanceSplittingFlows];*)
        NoDeadStarts = Join[outAuxEntryPairs, outAuxExitPairs, pairs];
        RuleBalanceGatheringFlows = Association[(j@@# -> Total[j@@@ FlowGathering[auxTriples][#]]) & /@ NoDeadStarts];(*Rule*)
        
        (*get equations for the exit currents at the entry vertices*)
        BalanceGatheringFlows = ((-j@@# + Total[j@@@ FlowGathering[auxTriples][#]]) & /@ NoDeadStarts);
        EqBalanceGatheringFlows = Simplify /@ (And @@ (# == 0 & /@ BalanceGatheringFlows));
        (*Print["nds: ",NoDeadStarts,"\n", RuleBalanceGatheringFlows, BalanceGatheringFlows, EqBalanceGatheringFlows];*)
        
        (*Incoming currents*)
        EqEntryIn = (j@@# == EntryDataAssociation[#]) & /@  inAuxEntryPairs;(*List of Equals*)
        
        (*Outgoing flows at entrances and incoming flows at exits are zero*)
        RuleEntryOut = Association[(j@@# -> 0) & /@ outAuxEntryPairs];
        RuleExitFlowsIn = Association[(j@@# -> 0) & /@ inAuxExitPairs];
		
		(*Exit values at exit vertices*)
        (*TODO: continue from here*)
        (*Print[ExitCosts,exitEdges];*)
        (*The value function on the auxiliary edges is constant and equal to the exit cost.*)
        RuleExitValues = AssociationThread[u@@@(Reverse/@outAuxExitPairs),Last /@ exitVerticesCosts];
        RuleExitValues = Join[RuleExitValues, AssociationThread[u@@@outAuxExitPairs, Last /@ exitVerticesCosts]];
        RuleEntryValues = AssociationThread[u@@@outAuxEntryPairs, u@@@inAuxEntryPairs];
        (*Print[RuleEntryValues];*)
        (*Print["new exit: ", RuleExitValues];*)
        EqValueAuxiliaryEdges = And @@ ((u@@# == u@@Reverse[#]) & /@ Join[inAuxEntryPairs, outAuxExitPairs]);
        (*Print[EqValueAuxiliaryEdges];*)
        (*use ToRules to get the rules*)
        (*Print[TransitionsAt[auxiliaryGraph, #] & /@ verticesList];*)
        IneqSwitchingByVertex = IneqSwitch[u, SwitchingCosts] @@@ TransitionsAt[auxiliaryGraph, #] & /@ verticesList;
        (*Print[IneqSwitchingByVertex];*)
        IneqSwitchingByVertex = And @@@ IneqSwitchingByVertex;
        (*Print[IneqSwitchingByVertex];*)
        (*IneqSwitchingByVertex = Select[#,FreeQ[Infinity]]&/@IneqSwitchingByVertex;
        Print[IneqSwitchingByVertex];*)
        (*Print[SwitchingCosts];*)
        AltOptCond = And @@ AltSwitch[j, u, SwitchingCosts] @@@ auxTriples;(*Or*)
        (*Print[AltOptCond];*)
        (*Print[SignedFlows];*)
        Nlhs = Flatten[u@@# - u@@Reverse@# + SignedFlows[#] & /@ halfPairs];
        (*Print[Nlhs];*)
        Nrhs = Flatten[SignedFlows[#] - Sign[SignedFlows[#]] Cost[SignedFlows[#], #] & /@ halfPairs];
        (*Print[Nrhs];*)
       	(*stuff to solve the general case faster*)
       	(*TODO: maybe change the notation to include the vertices*)
        costpluscurrents = Table[Symbol["cpc" <> ToString[k]], {k, 1, Length@edgeList}];
        EqGeneral = And @@ (MapThread[Equal, {Nlhs, costpluscurrents}]);
        costpluscurrents = AssociationThread[costpluscurrents, Nrhs];
        (*Print[EqGeneral, costpluscurrents];*)
        (*stuff to solve the general case faster*)

        (*list of all module variables, except for ModuleVars*)
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
          OutRules, InRules, IneqSwitchingByVertex, AltOptCond, Nlhs, 
          CostArgs, Nrhs, costpluscurrents, EqGeneral};
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
          "RuleExitValues", "EqValueAuxiliaryEdges", "OutRules", "InRules",
           "IneqSwitchingByVertex", "AltOptCond", "Nlhs", "CostArgs", 
          "Nrhs", "costpluscurrents", "EqGeneral"};
        Join[Data, AssociationThread[ModuleVarsNames, ModuleVars]](******)
        (**)
    ];
NumberVectorQ::usage = 
"NumberVectorQ[j] returns True if the vetor j is numeric."
NumberVectorQ[j_] :=
    And @@ (NumberQ /@ j);

GetKirchhoffMatrix::usage = "GetKirchhoffMatrix[d2e] returns the \
entry current vector, Kirchhoff matrix,  (critical congestion) cost \
function, and the variables in the order corresponding to the Kirchhoff matrix."
(*TODO: remove jvars, jtvars from code.*)
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
    
MFGPreprocessing::usage =
"MFGPreprocessing[Eqs] returns the association Eqs with the preliminary solution \"InitRules\" and corresponding \'reduced\' \"NewSystem\"."
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
        
        Print[Eqs["auxiliaryGraph"]];
        (*First rules:entry currents*)
        InitRules = Association[Flatten[ToRules /@ EqEntryIn]];
        (*no exit at the entrances*)
        AssociateTo[InitRules, RuleEntryOut];
        (*no entrance at the exits*)
        AssociateTo[InitRules, RuleExitFlowsIn];
        (*value function:  auxiliary edges*)
        AssociateTo[InitRules, RuleEntryValues];
        (*currents gathered from transition currents*)
        AssociateTo[InitRules, RuleBalanceGatheringFlows];
        (*value function: exit costs*)
        AssociateTo[InitRules, RuleExitValues];
        (*Print["InitRules",InitRules];*)
        temp = PrintTemporary["Simplifying inequalities involving Switching Costs...",IneqSwitchingByVertex];
        {time, IneqSwitchingByVertex} = AbsoluteTiming[And @@ (Simplify /@ ( IneqSwitchingByVertex /. InitRules))];
        NotebookDelete[temp];
        Print["Switching costs simplified in ", time, " seconds."(*,"\n", IneqSwitchingByVertex*)];
        EqBalanceSplittingFlows = EqBalanceSplittingFlows /. InitRules;
        EqValueAuxiliaryEdges = EqValueAuxiliaryEdges /. InitRules;
        temp = PrintTemporary["Solving some balance equations: ", EqBalanceSplittingFlows && EqValueAuxiliaryEdges];
        {time, Rules} = AbsoluteTiming[First@Solve[EqBalanceSplittingFlows && EqValueAuxiliaryEdges] // Quiet];
        InitRules = Join[InitRules /. Rules, Association@Rules];
        NotebookDelete[temp];
        Print["Balance equations solved in ", time, " seconds."];
        NewSystem = 
         MapThread[
          And, {{True, IneqJts && IneqJs, AltFlows && AltTransitionFlows && AltOptCond}, 
           Sys2Triple[IneqSwitchingByVertex]}
         ];
        temp = PrintTemporary["TripleClean will work on\n", NewSystem,"\nwith: ",InitRules];
        {time,{NewSystem, InitRules}} = AbsoluteTiming@TripleClean[{NewSystem, InitRules}];
        NotebookDelete[temp];
        Print["TripleClean took ", time, " seconds."];
        (*NewSystem[[2]] = DeleteDuplicates[NewSystem[[2]]];
        NewSystem[[3]] = DeleteDuplicates[NewSystem[[3]]];*)
        NewSystem[[1]] = EqGeneral;
        temp = PrintTemporary["TripleClean again (with some more equalities)..."(*,"\n", EqGeneral*)];
        {time,{NewSystem, InitRules}} = AbsoluteTiming@TripleClean[{NewSystem, InitRules}];
        NotebookDelete[temp];
        Print["TripleClean again (with some more equalities) took ", time, " seconds."];
        
        (*Print["Preprocessing: NewSystem after second TripleClean: ", NewSystem];*)
(*        NewSystem[[2]] = DeleteDuplicates[Simplify/@NewSystem[[2]]];
        NewSystem[[3]] = DeleteDuplicates[Simplify/@NewSystem[[3]]];
        Print["Preprocessing: NewSystem after second TripleClean deleteduplicates: ", NewSystem];*)
        ModuleVarsNames = {"InitRules", "NewSystem"};
        ModulesVars = {InitRules, NewSystem};
        Join[Eqs, AssociationThread[ModuleVarsNames, ModulesVars]]
    ];

CriticalCongestionSolver::usage = 
  "CriticalCongestionSolver[Eqs] returns Eqs with an association \"AssoCritical\" with rules to 
the solution to the critical congestion case";
CriticalCongestionSolver[$Failed] :=
    $Failed

CriticalCongestionSolver[Eqs_] :=
    Module[ {PreEqs, js, AssoCritical, time, temp},
    	If[KeyExistsQ[Eqs,"InitRules"],
    		PreEqs = Eqs,
    		temp = PrintTemporary["Preprocessing..."];
        	{time, PreEqs} = AbsoluteTiming@MFGPreprocessing[Eqs];
        	NotebookDelete[temp];
        	Print["Preprocessing took ", time, " seconds to terminate."];
        ];
        js = Lookup[PreEqs, "js", $Failed];
        AssoCritical = MFGSystemSolver[PreEqs][AssociationThread[js, 0 js]];
        Join[PreEqs, Association["AssoCritical"-> AssoCritical]]
    ];

MFGSystemSolver::usage = 
"MFGSystemSolver[Eqs][edgeEquations] returns the
association with rules to the solution";
MFGSystemSolver[Eqs_][approxJs_] :=
    Module[ {NewSystem, InitRules, pickOne, vars, System, Ncpc,
        costpluscurrents, us, js, jts, jjtsR, usR, time, temp,
        systemByTransition, ineqsByTransition, altsByTransition},
        us = Lookup[Eqs, "us", $Failed];
        js = Lookup[Eqs, "js", $Failed];
        jts = Lookup[Eqs, "jts", $Failed];
        InitRules = Lookup[Eqs, "InitRules", $Failed];
        NewSystem = Lookup[Eqs, "NewSystem", $Failed];
        costpluscurrents = Lookup[Eqs, "costpluscurrents", $Failed];
        {time, Ncpc} = AbsoluteTiming[RoundValues @ (Expand/@(costpluscurrents /. approxJs))];
        Print["MFGSS: Calculated the cost plus currents for the flow in ", time, " seconds."(*,"\n", Ncpc*)];
        InitRules = Expand/@(InitRules /. Ncpc);
        (*Print[NewSystem];*)
        (*Print["I know that for this example the inequalities are returning False:\n", NewSystem[[2]]];*)
        NewSystem = NewSystem /. Ncpc;
        (*Print["MFGSS: Numeric right-hand side:\n", NewSystem];*)
        
        If[And@@((#===False)&/@NewSystem),
        	Print["MFGSS: There is no solution!"];
            Return[Null,Module]
        ];
        (*Retrieve some equalities from the inequalites: group by transition flow.*)
        temp = PrintTemporary["MFGS: Selecting inequalities by transition flow..."];
        {time, ineqsByTransition} = AbsoluteTiming[Select[NewSystem[[2]], Function[exp, !FreeQ[#][exp]]]&/@jts];
        NotebookDelete[temp];
        Print["MFGS: Selecting inequalities by transition flow took ", time, " seconds. ", Length[ineqsByTransition]];
        temp = PrintTemporary["MFGSS: Simplifying inequalities by transition flow..."];
        {time, ineqsByTransition}= AbsoluteTiming[RemoveDuplicates[Simplify/@ineqsByTransition]];
        NotebookDelete[temp];
        Print["MFGSS: Simplifying inequalities by transition flow took ", time, " seconds. "(*, DeleteDuplicates[Select[#,(Head[#]===Equal)&]&/@DeleteCases[ineqsByTransition,True]]*)(*,"\n", ineqsByTransition*)];
		
		NewSystem[[2]]=(*Simplify@*)(And@@ineqsByTransition);
		NewSystem = Sys2Triple[And@@NewSystem];
		(*Print["MFGSS: TripleClean"];*)
		{NewSystem,InitRules} = TripleClean[{NewSystem,InitRules}];
		(*Print["TripleClean finished.","\nHere is the system:\n",NewSystem,"\nand here are the rules\n",InitRules];*)



		(*Return[{NewSystem,InitRules},Module];*)



		(*Print[NewSystem[[2]], "\n",AbsoluteTiming[Simplify@NewSystem[[2]]],"\n",AbsoluteTiming[Simplify/@NewSystem[[2]]]];*)
		
        (*Print["There are no equations: ", NewSystem[[1]]];
        temp = PrintTemporary["MFGS: Selecting inequalities by transition flow..."];
        {time, ineqsByTransition} = AbsoluteTiming[Select[NewSystem[[2]], Function[exp, !FreeQ[#][exp]]]&/@jts];
        NotebookDelete[temp];
        Print["MFGS: Selecting inequalities by transition flow took ", time, " seconds."];
        temp = PrintTemporary["MFGSS: Simplifying inequalities by transition flow..."];
        {time, ineqsByTransition}= AbsoluteTiming[Simplify/@ineqsByTransition];
        NotebookDelete[temp];
        Print["MFGSS: Simplifying inequalities by transition flow took ", time, " seconds.","\n", ineqsByTransition];*)
		
		
		
		(*temp = PrintTemporary["MFGS: Selecting alternatives by transition flows..."];
        {time, altsByTransition} = AbsoluteTiming[Select[NewSystem[[3]], Function[exp, !FreeQ[#][exp]]]&/@jts];
        NotebookDelete[temp];
        Print["MFGS: Selecting alternatives by transition flow took ", time, " seconds.","\n",altsByTransition];
        temp = PrintTemporary["MFGSS: Simplifying by transition flow..."];
        systemByTransition = MapThread[And,{ineqsByTransition, altsByTransition}];
        {time, systemByTransition}= AbsoluteTiming[Simplify/@systemByTransition];
        NotebookDelete[temp];
        Print["MFGSS: Simplifying by transition flow took ", time, " seconds.","\n", systemByTransition];
        NewSystem = And@@systemByTransition;  
        NewSystem = Sys2Triple[And@@NewSystem];
        {NewSystem,InitRules} = TripleClean[{NewSystem,InitRules}];*)
        (*Print["MFGSS: Simplifying\n", NewSystem[[2]]];
        NewSystem[[2]] = Simplify@NewSystem[[2]];*)
        
        (*Let's see if this is just slowing us down...
        If[NewSystem[[3]]=!=True, 
        	NewSystem[[3]] = SortBy[Simplify`SimplifyCount][Simplify/@NewSystem[[3]] ]
        ];*)
        
        
        {NewSystem, InitRules} = FinalStep[{NewSystem, InitRules}];
        
        
        System = BooleanConvert[NewSystem[[3]]]&&NewSystem[[2]];
        Print["BooleanConvert done. Simplifying..."];
        System = Simplify[And@@NewSystem];
        Print["Simplifying done. TripleClean..."];
        {System, InitRules} = TripleClean[{Sys2Triple[System], InitRules}];
        (*vars = Join[jts, js, us];
        vars = Intersection[Flatten[List@@@((Join[us, js] /. InitRules)/.Times->Plus)], vars];
        If[vars ==={},
        	Print["MFGSS: System is ", And@@NewSystem];
        	Return[InitRules, Module],
        	Print["MFGSS: These variables still need to be determined: ",vars];
        	(*Print["The system is:\n",NewSystem, "\nAnd the rules are:\n",InitRules];*)
        	(*Return[Null, Module]*)
        ];*)
        System = And@@System;
        Which[System === False, 
            Print["MFGSS: There is no solution!"];
            Return[Null,Module], 
            (*NewSystem[[1]]&&NewSystem[[3]] === True &&*) 
            System =!= True,
            Print["MFGSS: (Possibly) Multiple solutions:\n", System];
            jjtsR = Join[js, jts];
            jjtsR = Intersection[Flatten[List@@@((Join[us, js] /. InitRules)/.Times->Plus)], jjtsR];
        	usR = Intersection[Flatten[List@@@((Join[us, js] /. InitRules)/.Times->Plus)], us];
            
            (*usR = Select[us, Not[FreeQ[System, #]] &];
            jjtsR = Select[Join[js, jts], Not[FreeQ[System, #]] &];*)
            vars = Join[usR, jjtsR];
            
            (*Have to pick one so that all the currents have numerical values*)
            pickOne = FindInstance[System && And @@ ((# > 0 )& /@ jjtsR), vars, Reals];
            If[pickOne === {}, 
            	pickOne = FindInstance[System && And @@ ((# >= 0 )& /@ jjtsR), vars, Reals]
            ];
           	pickOne = Association @ First @ pickOne;
            Print["MFGSS: Picked one solution: ", pickOne];
            InitRules = Expand /@ Join[InitRules /. pickOne, pickOne],
            True, 
         	Print["MFGSS: System is ", System](*Reducing in the reals should be a good way of simplifying things here.*)
         ];
        InitRules = Join[KeyTake[InitRules, us], KeyTake[InitRules, js], KeyTake[InitRules, jts]];
        InitRules
    ];


FinalStep::usage = 
"FinalStep[{EE,NN,OR}, rules] takes a grouped system and some Association of rules (a partial solution). It returns the result of applying ZAnd ";

FinalStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, rules_}] :=
    {{EE, NN, OR}, rules}

FinalStep[{{EE_, NN_, OO_}, rules_}] :=
    Module[ {NewSystem, newrules, sorted = True, time, temp, NewSystemBC},
        {NewSystem, newrules} = TripleClean[{{EE, NN, OO}, rules}];
        If[NewSystem[[3]] === True,
        	(*only inequalities*)
            Return[{NewSystem, newrules}, Module],
            Print["Final: ", Length/@NewSystem];
            (*temp = PrintTemporary["Final: RemoveDuplicated alternatives..."] ;*)
            (*{time, sorted} = AbsoluteTiming@RemoveDuplicates[NewSystem[[3]]];*)
            sorted = RemoveDuplicates[NewSystem[[3]]];
            (*Print["Final: sorted: ",sorted];*)
            (*NotebookDelete[temp];
            Print["Final: RemoveDuplicated alternatives took ", time, " seconds. ", Length[sorted]];*)
        ];
        (*Print["Final:\n", NewSystem, "\nRules:\n", newrules];*)
        temp = PrintTemporary["Final: Iterative DNF convertion on " , Length[sorted]," disjunctions..."(*,  Append[sorted][Most[NewSystem]]*)];
        {time, NewSystem} = AbsoluteTiming[ZAnd[And @@ Most[NewSystem], sorted]];
        (*Print[Length[NewSystem], " ", Length/@NewSystem];*)
        NotebookDelete[temp];
        Print["Final: Iterative DNF convertion on " , Length[sorted]," disjunctions took ", time, " seconds to terminate."(*,"\nGiving\n", NewSystem*)];
        If[Head[NewSystem] === Or, 
	        (*temp =*) Print(*Temporary*)["Final: Reducing (each alternative)... ", Length[NewSystem] ," of them.","\n", Head[NewSystem],"\n", NewSystem] ;
        	(*Print["Final: Head is OR: ", NewSystem];*)
        	{time, NewSystem} = AbsoluteTiming[Reduce[#,Reals]&/@NewSystem];
        	NewSystem = Simplify[NewSystem],
        	(*Print["Head is not OR. Should we BooleanConvert anyway?","\n", NewSystem];*)
        	{time,NewSystem} = AbsoluteTiming@Reduce[NewSystem, Reals]
        ];
        (*NotebookDelete[temp];*)
        Print["Final: Reducing (each alternative), returned ", Length@NewSystem," of them, ","took ", time, " seconds to terminate."];
        (*NewSystem = Simplify@NewSystem;*)
        NewSystem = BooleanConvert@NewSystem;
        (*Print["Final: ", NewSystem];*)
        NewSystem = Sys2Triple[NewSystem];
        {NewSystem, newrules} = TripleClean[{NewSystem, newrules}];
        Print["Now: ", TimeObject[Now], "The new rules are: ", newrules, Length/@NewSystem];
        {NewSystem, newrules}
    ];


(*FinalClean::usage =
"FinalClean[{{EE,NN,OR},rules}] performs FinalStep, Reduce, and TripleClean on {{EE,NN,OR},rules}.
The result should always be {{True, True, NewOr}, NewRules}"
FinalClean[{{EE_, NN_, OR_}, rules_}] := 
With[
	{NewSystemRules = FinalStep[{{EE, NN, OR}, rules}]},
	(*Print["FinalClean... ", NewSystemRules[[1]]];*)
	TripleClean[{(Simplify/@NewSystemRules[[1]]),NewSystemRules[[2]]}]
];*)

Sys2Triple::usage =
"Sys2Triple[sys] returns a triple with equalities, inequalites, and alternatives from sys, respectively.
The input, sys, should be a system of equations, inequalities and (simple) alternatives."
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
          (*This works because of the structure of our system. *)
          NN = And @@ Lookup[groups, False, True];
          {EE, NN, OR}
      ], 
     Head[system] === Equal, 
      {system, True, True}, 
     Head[system] === Or, 
         {True, True, system}, 
     True, 
         {True, system, True}];

TripleStep::usage =
"TripleStep[{{EE,NN,OR},Rules}] returns {{NewEE, NewNN, NewOR}, NewRules}, where NewRules contain the solutions to all the equalities found in the system 
after replacing Rules in {EE,NN,OR}."

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
	(*Print["TripleStep 1: ", EE];*)
	If[ EE =!= True && EE =!= False,
    	newrules = First@Solve[EE] // Quiet
    ];
    newrules = Expand /@ Join[rules /. newrules, Association@newrules];
    (*Print["TripleStep 1.3: ", NNs /. newrules];*)
    NN = Expand /@ (NNs /. newrules);
    (*Print["TripleStep 1.5: ", NN /. rules];*)
    OR = Expand /@ (ORs /. newrules);
    (*Print["TripleStep 2: ", OR];*)
(*    NN = Simplify[NN];
     Print["TripleStep 3: ", NN];*)
    {NNE, NN, NNO} = Sys2Triple[NN];
    {ORE, ORN, OR} = Sys2Triple[OR];
    EE = NNE && ORE;
    (*NNO and ORN are empty because no OR is generated from inequalities and vice-versa.*)
    {{EE, NN, OR}, newrules}
    ];


(*TripleStep[{{EEs_, NNs_, ORs_}, rules_Association}] :=
    Module[ {EE = EEs /. rules, NN = Simplify /@ (NNs /. rules), 
      OR = Simplify /@ (ORs /. rules), NNE, NNO, ORE, ORN,
      newrules = {}},
      Print["TripleStep 1: ",{EE, NN,OR}];
        NN = Simplify[NN];
      Print["TripleStep 2: ", NN];
        {NNE, NN, NNO} = Sys2Triple[NN];
        {ORE, ORN, OR} = Sys2Triple[OR];
        EE = EE && NNE && ORE;
        If[ EE =!= True && EE =!= False,
            newrules = First@Solve[EE] // Quiet
        ];
        newrules = Expand /@ Join[rules /. newrules, Association@newrules];
        {{EE, NN, OR}, newrules}
    ];*)

TripleClean::usage =
"TripleClean[{{EE,NN,OR},Rules}] composes TripleStep until it reaches a fixed point, that is, {{True,NewNN,NewOR},NewRules} such that 
replacement of NewRules in NewNN and NewOR do not produce equalities."
TripleClean[{{EE_, NN_, OR_}, rules_}] := FixedPoint[TripleStep, {{EE, NN, OR}, rules}];

ZAnd::usage =
"
ZAnd[xp,xps] returns a system which is equivalent to xp&&xps in disjunctive normal form. 
ZAnd[xp, And[fst, scd] 
ZAnd[xp, eq] returns xp with the solution of eq replaced in it together with eq.
ZAnd[xp, Or[fst,scd]] returns ZAnd[xp, fst]||ZAnd[xp, scd]";

ZAnd[_, False] :=
    False

ZAnd[False, _] :=
    False

ZAnd[xp_, True] :=
    xp

ZAnd[xp_, eq_Equal] := ReZAnd[xp, True, eq]

ZAnd[xp_, And[fst_,rst_]] :=
    If[ Head[fst] === Or,
        RemoveDuplicates@(ReZAnd[xp, rst] /@ fst),
        ReZAnd[xp, rst, fst]
    ]

ZAnd[xp_, Or[fst_,scd_]] :=
	With[{rfst = Reduce[fst,Reals]},
    	RemoveDuplicates@(Or@@(ZAnd[xp, #] & /@ {rfst,scd}))
	];

ZAnd[xp_, leq_] := xp && leq

(*Operator form of ReZAnd*)
ReZAnd[xp_, rst_] :=
    ReZAnd[xp, rst, #] &

(*ReZAnd[xp_, rst_, fst_Equal] :=
        If[ Simplify[fst] === False,
            False,
            With[ {fsol = First@Solve@fst},
                ZAnd[Simplify[(xp /. fsol)] && fst, ReplaceSolution[rst, fsol]]
            ]
        ]*)

ReZAnd[xp_, rst_, fst_Equal] :=
Module[{newfst=Simplify@fst},
        If[ newfst === False,
            False,
            With[ {fsol = First@Solve@newfst},
                ZAnd[ReplaceSolution[xp, fsol] && fst, ReplaceSolution[rst, fsol]]
            ]
        ]
]


ReZAnd[xp_, rst_, fst_] :=
    ZAnd[xp && fst, rst]

Clear[ReplaceSolution];
ReplaceSolution::usage =
"ReplaceSolution[xp,sol] substitutes the Rule, solution, on the expression xp. After that, it simplifies the first equation if the Head of the expression is And. Otherwise, it simplifies the whole expression.";
ReplaceSolution[rst_?BooleanQ, sol_] :=
    rst


ReplaceSolution[rst_, sol_] :=
    With[ {newrst = rst /. sol},
    	If[ Head[newrst] === And,
            And[Simplify@First@newrst, Rest@newrst],
            Simplify[newrst]
        ]
    ]
    

(*ReplaceSolution[rst_And, sol_] :=
    Module[{groups=GroupBy[List@@rst, (!FreeQ[And@@(First/@sol)][#])&], reduced},
    	(*Print["reducing "];*)
    	(*reduced = Reduce[(And@@Lookup[groups, True, {}]/.sol),Reals];*)
    	(*Print[Lookup[groups, True, {}]/.sol];*)
    	reduced = And@@RemoveDuplicates[Lookup[groups, True, {}]/.sol];
    	(*reduced = Simplify[reduced];
    	Print["reduced: ", reduced];*)
    	reduced&&And@@Lookup[groups, False, {}]
    	])

ReplaceSolution[rst_, sol_] := 
Reduce[rst/.sol, Reals]
*)
SortOp = SortBy[Simplify`SimplifyCount]

RemoveDuplicates::usage =
"RemoveDuplicates[xp] sorts and then DeleteDuplicates. 
We need to sort because DeleteDuplicates only deletes identical expressions.
For example (A&&B)||(B&&A) becomes (A&&B) only after sorting.
"
RemoveDuplicates[xp_And] :=
    DeleteDuplicates[SortOp[xp]];

RemoveDuplicates[xp_Or] :=
    DeleteDuplicates[SortOp[xp]];

RemoveDuplicates[xp_] :=
    xp