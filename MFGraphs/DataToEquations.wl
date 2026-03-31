(*Wolfram Language package*)

(* --- Public API declarations --- *)

ConsistentSwitchingCosts::usage =
"ConsistentSwitchingCosts[switchingcosts][{a,b,c}->S]
returns True if S, the cost of switching from the edge ab to cb, is smaller than any other combination,
such as, ab to bd and then from db to bc.
Returns the condition for this switching cost to satisfy the triangle inequality when S, and the other
switching costs too, does not have a numerical value.";

IsSwitchingCostConsistent::usage =
"IsSwitchingCostConsistent[List of switching costs] is True if all switching costs satisfy the triangle inequality. If some switching costs are symbolic, then it returns the consistency conditions."

AltFlowOp::usage =
"AltFlowOp[j][list] returns the alternative: j@@list ==0 || j@@Reverse@list ==0.";

FlowSplitting::usage =
"FlowSplitting[AT][UndirectedEdge[a, b]] returns the splitting that start with {a,b}.";

FlowGathering::usage=
"FlowGathering[auxTriples_List][x_] returns the triples that end with x.";

IneqSwitch::usage =
"IneqSwitch[u, switchingCosts][{v,e1,e2}] returns the optimality condition at the vertex v related to switching from e1 to e2. Namely,
u[{v, e1}] <= u[{v, e2}] + switchingCosts[{v, e1, e2}]"

AltSwitch::usage =
"AltSwitch[jt,u,switchingCosts][{v,e1,e2}] returns the complementarity condition:
(jt[{v, e1, e2}] == 0) || (u[{v, e2}] - u[{v, e1}] + switchingCosts[{v, e1, e2}] == 0)"

DataToEquations::usage = "DataToEquations[Data] returns the equations, \
inequalities, and alternatives associated to the Data. "

NumberVectorQ::usage =
"NumberVectorQ[j] returns True if the vector j is numeric."

GetKirchhoffMatrix::usage = "GetKirchhoffMatrix[d2e] returns the \
entry current vector, Kirchhoff matrix, (critical congestion) cost \
function placeholder, and the variables in the order corresponding to the Kirchhoff matrix. \
The third slot is retained for backward compatibility and should not be used by new code."

GetKirchhoffLinearSystem::usage =
"GetKirchhoffLinearSystem[d2e] returns the entry current vector, Kirchhoff matrix, \
and the variables in the order corresponding to the Kirchhoff matrix.";

MFGPreprocessing::usage =
"MFGPreprocessing[Eqs] returns the association Eqs with the preliminary solution \"InitRules\" and corresponding 'reduced' \"NewSystem\"."

CriticalCongestionSolver::usage =
  "CriticalCongestionSolver[Eqs] returns Eqs with an association \"AssoCritical\" with rules to
the solution to the critical congestion case. Option: \"ReturnShape\" (default \"Legacy\"); \
use \"Standard\" to add normalized solver-result keys.";

MFGSystemSolver::usage =
"MFGSystemSolver[Eqs][edgeEquations] returns the
association with rules to the solution";

DNFSolveStep::usage =
"DNFSolveStep[{EE,NN,OR}, rules] takes a grouped system and some Association of rules (a partial solution). It returns the result of applying DNFReduce.";

SystemToTriple::usage =
"SystemToTriple[sys] returns a triple {equalities, inequalities, alternatives} from sys.
The input should be a system of equations, inequalities and (simple) alternatives."

TripleStep::usage =
"TripleStep[{{EE,NN,OR},Rules}] returns {{NewEE, NewNN, NewOR}, NewRules}, where NewRules contain the solutions to all the equalities found in the system after replacing Rules in {EE,NN,OR}."

TripleClean::usage =
"TripleClean[{{EE,NN,OR},Rules}] composes TripleStep until it reaches a fixed point, that is, {{True,NewNN,NewOR},NewRules} such that replacement of NewRules in NewNN and NewOR do not produce equalities."

Data2Equations::usage = "Data2Equations is a backward-compatibility alias for DataToEquations.";
FinalStep::usage = "FinalStep is a backward-compatibility alias for DNFSolveStep.";

DataToEquations::switchingcosts =
"Switching costs are inconsistent.";

MFGSystemSolver::nosolution =
"There is no feasible symbolic solution for the current system.";

Begin["`Private`"];

(* --- Switching cost consistency --- *)

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

IsSwitchingCostConsistent[switchingCosts_] :=
    And @@ Simplify[ConsistentSwitchingCosts[switchingCosts] /@ switchingCosts]

(* --- Graph helper functions --- *)

TransitionsAt[G_, k_] :=
    Insert[k,2] /@ Permutations[AdjacencyList[G, k], {2}]

AltFlowOp[j_][list_]:=
	j@@list==0||j@@Reverse@list==0;

FlowSplitting[auxTriples_List][x_] := Select[auxTriples, MatchQ[#,{Sequence@@x,__}]&] ;

FlowGathering[auxTriples_List][x_] := Select[auxTriples, MatchQ[#,{__,Sequence@@x}]&]

IneqSwitch[u_,Switching_Association][r_,i_,w_] := u[r,i] <= u[w,i]+Switching[{r,i,w}];

AltSwitch[j_, u_, Switching_][r_, i_, w_] :=
    (j[r,i,w] == 0) ||
    (u[r,i] == u[w,i] + Switching[{r,i,w}]);

(* --- DataToEquations: main converter --- *)

DataToEquations[Data_Association] :=
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
            Message[DataToEquations::switchingcosts];
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

GetKirchhoffLinearSystem[Eqs_] :=
    Module[ {Kirchhoff, EqEntryIn = Lookup[Eqs, "EqEntryIn", True],
      BalanceGatheringFlows =
       Lookup[Eqs, "BalanceGatheringFlows", {}],
      BalanceSplittingFlows =
       Lookup[Eqs, "BalanceSplittingFlows", {}],
      RuleExitFlowsIn = Lookup[Eqs, "RuleExitFlowsIn", {}],
      RuleEntryOut = Lookup[Eqs, "RuleEntryOut", {}], BM, KM, vars,
      cost, CCost},
        Kirchhoff = Join[EqEntryIn, (# == 0 & /@ (BalanceGatheringFlows + BalanceSplittingFlows))];
        Kirchhoff = Kirchhoff /. Join[RuleExitFlowsIn, RuleEntryOut];
        (* Extract all flow variables from the Kirchhoff system *)
        vars = Select[Variables[Kirchhoff /. Equal -> List], MatchQ[#, j[_,_,_] | j[_,_]] &];
        If[Length[vars] === 0,
          (* Degenerate case: no flow variables *)
          {0, 0, {}},
          {BM, KM} = CoefficientArrays[Kirchhoff, vars];
          {-BM, KM, vars}
        ]
    ];

GetKirchhoffMatrix[Eqs_] :=
    Module[{B, K, vars, cost, CCost},
        {B, K, vars} = GetKirchhoffLinearSystem[Eqs];
        If[Length[vars] === 0,
            Return[{B, K, <||>, vars}, Module]
        ];
        (* The third slot is a legacy zero-cost placeholder retained only for compatibility. *)
        cost = AssociationThread[vars, (0 &) /@ vars];
        CCost = cost /@ vars /. MapThread[Rule, {vars, #}] &;
        {B, K, CCost, vars}
    ];

(* --- MFGPreprocessing --- *)

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
        With[{items = IneqSwitchingByVertex /. InitRules},
          {time, IneqSwitchingByVertex} = AbsoluteTiming[
            And @@ If[Length[items] >= $MFGraphsParallelThreshold,
              (If[$KernelCount === 0, LaunchKernels[]]; ParallelMap[Simplify, items]),
              Simplify /@ items
            ]
          ]
        ];
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
           SystemToTriple[IneqSwitchingByVertex]}
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

Options[CriticalCongestionSolver] = {"ReturnShape" -> "Legacy"};

CriticalCongestionSolver[$Failed, ___] := $Failed

CriticalCongestionSolver[Eqs_, OptionsPattern[]] :=
    Module[ {PreEqs, js, AssoCritical, time, temp, status, returnShape,
        resultKind, message, solution},
        returnShape = OptionValue["ReturnShape"];
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
        (* Feasibility check: any flow variable with negative numeric value *)
        status = If[AssoCritical === Null, "Infeasible",
            Module[{flowKeys, flowVals},
                flowKeys = Select[Keys[AssoCritical], MatchQ[#, _j] &];
                flowVals = Lookup[AssoCritical, flowKeys];
                If[flowVals === {} || Min[Select[flowVals, NumericQ]] < 0, "Infeasible", "Feasible"]
            ]
        ];
        If[returnShape === "Standard",
            resultKind = If[AssoCritical === Null, "Failure", "Success"];
            message = If[AssoCritical === Null, "NoSolution", None];
            solution = If[AssociationQ[AssoCritical], AssoCritical, Missing["NotAvailable"]];
            Join[
                PreEqs,
                MakeSolverResult[
                    "CriticalCongestion",
                    resultKind,
                    status,
                    message,
                    solution,
                    <|"AssoCritical" -> AssoCritical, "Status" -> status|>
                ]
            ],
            Join[PreEqs, <|"AssoCritical" -> AssoCritical, "Status" -> status|>]
        ]
    ];

(* --- MFGSystemSolver --- *)

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
        	Message[MFGSystemSolver::nosolution];
            Return[Null,Module]
        ];

        (* Retrieve some equalities from the inequalities: group by transition flow *)
        temp = MFGPrintTemporary["MFGSS: Selecting inequalities by transition flow..."];
        If[NewSystem[[2]] === True,
            ineqsByTransition = ConstantArray[True, Length[jts]];
            time = 0.,
            {time, ineqsByTransition} = AbsoluteTiming[
              If[Length[jts] >= $MFGraphsParallelThreshold,
                (If[$KernelCount === 0, LaunchKernels[]];
                 With[{ineqs = NewSystem[[2]]},
                   ParallelMap[Function[jt, Select[ineqs, !FreeQ[jt][#]&]], jts]
                 ]),
                Select[NewSystem[[2]], Function[exp, !FreeQ[#][exp]]]&/@jts
              ]
            ]
        ];
        NotebookDelete[temp];
        MFGPrint["MFGSS: Selecting inequalities by transition flow took ", time, " seconds. ", Length[ineqsByTransition]];
        temp = MFGPrintTemporary["MFGSS: Simplifying inequalities by transition flow..."];
        {time, ineqsByTransition} = AbsoluteTiming[
          DeduplicateByComplexity[
            If[Length[ineqsByTransition] >= $MFGraphsParallelThreshold,
              (If[$KernelCount === 0, LaunchKernels[]];
               ParallelMap[Simplify, ineqsByTransition]),
              Simplify /@ ineqsByTransition
            ]
          ]
        ];
        NotebookDelete[temp];
        MFGPrint["MFGSS: Simplifying inequalities by transition flow took ", time, " seconds. "];

		NewSystem[[2]]=(And@@ineqsByTransition);
		NewSystem = SystemToTriple[And@@NewSystem];
		{NewSystem,InitRules} = TripleClean[{NewSystem,InitRules}];

		{NewSystem, InitRules} = DNFSolveStep[{NewSystem, InitRules}];

        System = BooleanConvert[NewSystem[[3]]]&&NewSystem[[2]];
        MFGPrint["BooleanConvert done. Simplifying..."];
        System = Simplify[And@@NewSystem];
        MFGPrint["Simplifying done. TripleClean..."];
        {System, InitRules} = TripleClean[{SystemToTriple[System], InitRules}];
        System = And@@System;
        Which[System === False,
            Message[MFGSystemSolver::nosolution];
            Return[Null,Module],
            System =!= True,
            MFGPrint["MFGSS: (Possibly) Multiple solutions:\n", System];
            (* Extract all unsolved variables from System *)
            vars = Select[Variables[System], MatchQ[#, j[_,_,_] | j[_,_] | u[_,_] | u[_,_,_]] &];
            vars = Complement[vars, Keys[InitRules]];
            jjtsR = Select[vars, MatchQ[#, j[_,_,_] | j[_,_]] &];
            usR = Select[vars, MatchQ[#, u[_,_,_] | u[_,_]] &];

            (* Pick one solution so that all the currents have numerical values *)
            If[Length[vars] > 0,
              (* Attempt to find a solution *)
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
        (* Resolve transitive rule chains, e.g. j[2,4] -> j[1,2,4]+j[3,2,4]
           with j[1,2,4] -> 50 becomes j[2,4] -> 50 *)
        InitRules = Expand /@ FixedPoint[Function[r, ReplaceAll[r] /@ r], InitRules, 10];
        InitRules = Join[KeyTake[InitRules, us], KeyTake[InitRules, js], KeyTake[InitRules, jts]];
        InitRules
    ];

(* --- DNFSolveStep --- *)

DNFSolveStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, rules_}] :=
    {{EE, NN, OR}, rules}

DNFSolveStep[{{EE_, NN_, OO_}, rules_}] :=
    Module[ {NewSystem, newrules, sorted = True, time, temp},
        {NewSystem, newrules} = TripleClean[{{EE, NN, OO}, rules}];
        If[NewSystem[[3]] === True,
            Return[{NewSystem, newrules}, Module],
            MFGPrint["Final: ", Length/@NewSystem];
            sorted = DeduplicateByComplexity[NewSystem[[3]]];
        ];
        temp = MFGPrintTemporary["Final: Iterative DNF conversion on " , Length[sorted]," disjunctions..."];
        {time, NewSystem} = AbsoluteTiming[DNFReduce[And @@ Most[NewSystem], sorted]];
        NotebookDelete[temp];
        MFGPrint["Final: Iterative DNF conversion on " , Length[sorted]," disjunctions took ", time, " seconds to terminate."];
        If[Head[NewSystem] === Or,
	        MFGPrint["Final: Reducing (each alternative)... ", Length[NewSystem] ," of them."];
        	{time, NewSystem} = AbsoluteTiming[
        	  If[Length[NewSystem] >= $MFGraphsParallelThreshold,
        	    (If[$KernelCount === 0, LaunchKernels[]];
        	     ParallelMap[Function[branch, Reduce[branch, Reals]], List @@ NewSystem]),
        	    Reduce[#, Reals]& /@ NewSystem
        	  ]
        	];
        	NewSystem = Simplify[Or @@ NewSystem],
        	{time,NewSystem} = AbsoluteTiming@Reduce[NewSystem, Reals]
        ];
        MFGPrint["Final: Reducing (each alternative), returned ", Length@NewSystem," of them, ","took ", time, " seconds to terminate."];
        NewSystem = BooleanConvert@NewSystem;
        NewSystem = SystemToTriple[NewSystem];
        {NewSystem, newrules} = TripleClean[{NewSystem, newrules}];
        MFGPrint["Now: ", TimeObject[Now], " The new rules are: ", newrules, Length/@NewSystem];
        {NewSystem, newrules}
    ];

(* --- SystemToTriple: decompose into equalities, inequalities, alternatives --- *)

SystemToTriple[True] = Table[True, 3]

SystemToTriple[False] = Table[False, 3]

SystemToTriple[system_] :=
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

(* SafeFirstSolve: returns a list of Rules from CachedSolve, or {} on failure *)
SafeFirstSolve[eq_] := Module[{sol},
    sol = Quiet[CachedSolve[eq]];
    If[MatchQ[sol, {{__Rule}, ___}], First[sol], {}]
];

TripleStep[{{EEs_, NNs_, ORs_}, rules_List}] :=
    TripleStep[{{EEs, NNs, ORs}, Association@rules}]

TripleStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, rules_Association}] :=
    {{EE, NN, OR}, rules}

TripleStep[{{EE_, NN_?TrueQ, OR_?TrueQ}, rules_Association}] :=
    Module[ {newrules = {}},
        newrules = SafeFirstSolve[EE /. rules];
        newrules = Join[rules /. newrules, Association@newrules];
        {{True, NN, OR}, Expand /@ newrules}
    ];

TripleStep[{{True, NNs_, ORs_}, rules_Association}] :=
	{{True, NNs, ORs}, rules}

TripleStep[{{EEs_, NNs_, ORs_}, rules_Association}] :=
	Module[{EE = EEs /. rules, NN, OR, NNE, NNO, ORE, ORN, newrules ={}},
	If[ EE =!= True && EE =!= False,
    	newrules = SafeFirstSolve[EE]
    ];
    newrules = Expand /@ Join[rules /. newrules, Association@newrules];
    NN = Expand /@ (NNs /. newrules);
    OR = Expand /@ (ORs /. newrules);
    {NNE, NN, NNO} = SystemToTriple[NN];
    {ORE, ORN, OR} = SystemToTriple[OR];
    EE = NNE && ORE;
    {{EE, NN, OR}, newrules}
    ];

TripleClean[{{EE_, NN_, OR_}, rules_}] := FixedPoint[TripleStep, {{EE, NN, OR}, rules}];

(* --- Backward compatibility aliases --- *)
Data2Equations = DataToEquations;
FinalStep = DNFSolveStep;

End[];
