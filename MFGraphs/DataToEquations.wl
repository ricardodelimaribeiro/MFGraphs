(* Wolfram Language package *)
(*
   DataToEquations: Converts network topological data into Mean Field Game (MFG) equations.
   
   Theoretical Basis:
   - "First-order mean-field games on networks and Wardrop equilibrium" (Al Saleh et al., 2024)
   - "Hessian Riemannian Flow For Multi-Population Wardrop Equilibrium" (Bakaryan et al., 2025)

   This module constructs the symbolic system of equations (Hamilton-Jacobi, Fokker-Planck, 
   and Vertex-Transition constraints) required to find stationary Wardrop equilibria 
   on directed/undirected networks with congestion.
*)

(* --- Public API declarations --- *)

ConsistentSwitchingCosts::usage = "ConsistentSwitchingCosts[switchingcosts][{a,b,c}->S]
returns True if S, the cost of switching from the edge ab to cb, is smaller than any other combination,
such as, ab to bd and then from db to bc.
Returns the condition for this switching cost to satisfy the triangle inequality when S, and the other
switching costs too, does not have a numerical value.";

IsSwitchingCostConsistent::usage = "IsSwitchingCostConsistent[List of switching costs] is True if all switching costs satisfy the triangle inequality. If some switching costs are symbolic, then it returns the consistency conditions."

AltFlowOp::usage = "AltFlowOp[j][list] returns the alternative: j@@list ==0 || j@@Reverse@list ==0.";

FlowSplitting::usage = "FlowSplitting[AT][UndirectedEdge[a, b]] returns the splitting that start with {a,b}.";

FlowGathering::usage = "FlowGathering[auxTriples_List][x_] returns the triples that end with x.";

IneqSwitch::usage = "IneqSwitch[u, switchingCosts][v, e1, e2] returns the optimality condition at the vertex v related to switching from e1 to e2. Namely,
u[v, e1] <= u[v, e2] + switchingCosts[{v, e1, e2}]"

AltSwitch::usage = "AltSwitch[j, u, switchingCosts][v, e1, e2] returns the complementarity condition:
(j[v, e1, e2] == 0) || (u[v, e1] == u[e2, e1] + switchingCosts[{v, e1, e2}])"

DataToEquations::usage = "DataToEquations[Data] returns the equations, inequalities, and alternatives associated to the Data. "

NumberVectorQ::usage = "NumberVectorQ[j] returns True if the vector j is numeric."

GetKirchhoffMatrix::usage = "GetKirchhoffMatrix[d2e] returns the entry current vector, Kirchhoff matrix, (critical congestion) cost function placeholder, and the variables in the order corresponding to the Kirchhoff matrix. The third slot is retained for backward compatibility and should not be used by new code."

GetKirchhoffLinearSystem::usage = "GetKirchhoffLinearSystem[d2e] returns the entry current vector, Kirchhoff matrix, and the variables in the order corresponding to the Kirchhoff matrix.";

MFGPreprocessing::usage = "MFGPreprocessing[Eqs] returns the association Eqs with the preliminary solution \"InitRules\" and corresponding 'reduced' \"NewSystem\"."

CriticalCongestionSolver::usage = "CriticalCongestionSolver[Eqs] returns Eqs with an association \"AssoCritical\" with rules to
the solution to the critical congestion case. It returns a standardized association
containing solver metadata, feasibility, comparison fields, and the solver-specific
payload key \"AssoCritical\".";

MFGSystemSolver::usage = "MFGSystemSolver[Eqs][edgeEquations] returns the
association with rules to the solution";

DNFSolveStep::usage = "DNFSolveStep[{EE,NN,OR}, rules] takes a grouped system and some Association of rules (a partial solution). It returns the result of applying DNFReduce.";

SystemToTriple::usage = "SystemToTriple[sys] returns a triple {equalities, inequalities, alternatives} from sys.
The input should be a system of equations, inequalities and (simple) alternatives."

TripleStep::usage = "TripleStep[{{EE,NN,OR},Rules}] returns {{NewEE, NewNN, NewOR}, NewRules}, where NewRules contain the solutions to all the equalities found in the system after replacing Rules in {EE,NN,OR}."

TripleClean::usage = "TripleClean[{{EE,NN,OR},Rules}] composes TripleStep until it reaches a fixed point, that is, {{True,NewNN,NewOR},NewRules} such that replacement of NewRules in NewNN and NewOR do not produce equalities."

Data2Equations::usage = "Data2Equations is a backward-compatibility alias for DataToEquations.";

FinalStep::usage = "FinalStep is a backward-compatibility alias for DNFSolveStep.";

DataToEquations::switchingcosts = "Switching costs are inconsistent.";

MFGSystemSolver::nosolution = "There is no feasible symbolic solution for the current system.";

Begin["`Private`"];

(* --- Switching cost consistency --- *)

ConsistentSwitchingCosts[sc_][{a_, b_, c_} -> S_] :=
    Module[{origin, bounds},
        origin = Cases[sc, HoldPattern[{a, b, _} -> _]];
        origin = DeleteCases[origin, {a, b, c} -> S];
        If[origin =!= {},
            bounds = ((S <= Last[#] + Association[sc][{Part[First[#],
                 3], b, c}])& /@ origin);
            And @@ bounds
            ,
            True
        ]
    ];

IsSwitchingCostConsistent[switchingCosts_] :=
    And @@ Simplify[ConsistentSwitchingCosts[switchingCosts] /@ switchingCosts
        ]

(* --- Graph helper functions --- *)

TransitionsAt[G_, k_] :=
    Insert[k, 2] /@ Permutations[AdjacencyList[G, k], {2}]

AltFlowOp[j_][list_] :=
    j @@ list == 0 || j @@ Reverse @ list == 0;

FlowSplitting[auxTriples_List][x_] :=
    Select[auxTriples, MatchQ[#, {Sequence @@ x, __}]&];

FlowGathering[auxTriples_List][x_] :=
    Select[auxTriples, MatchQ[#, {__, Sequence @@ x}]&]

IneqSwitch[u_, Switching_Association][r_, i_, w_] :=
    u[r, i] <= u[w, i] + Switching[{r, i, w}];

AltSwitch[j_, u_, Switching_][r_, i_, w_] :=
    (j[r, i, w] == 0) || (u[r, i] == u[w, i] + Switching[{r, i, w}]);

(* --- Graph-distance-based Or-resolution --- *)
(* Resolves flow complementarity conditions (j[a,b]==0 || j[b,a]==0)
   using shortest-path distances from exit vertices. The vertex farther
   from the exit sends flow toward the closer vertex, setting the
   backward flow to zero. This eliminates most Or-conditions before
   DNFReduce, reducing exponential branching to O(V+E) graph traversal. *)

ResolveOneOr[Or[lhs_, rhs_], distToExit_Association] :=
    Which[
        (* Edge flow: j[a,b]==0 || j[b,a]==0 *)
        MatchQ[{lhs, rhs}, {_j == 0, _j == 0}] &&
            Length[lhs[[1]]] == 2 && Length[rhs[[1]]] == 2 &&
            lhs[[1, 1]] === rhs[[1, 2]] && lhs[[1, 2]] === rhs[[1, 1]],
        With[{da = Lookup[distToExit, lhs[[1, 1]], Infinity],
              db = Lookup[distToExit, lhs[[1, 2]], Infinity]},
            Which[da > db, rhs, da < db, lhs, True, $Failed]
        ],
        (* Transition flow: j[a,b,c]==0 || j[c,b,a]==0 *)
        MatchQ[{lhs, rhs}, {_j == 0, _j == 0}] &&
            Length[lhs[[1]]] == 3 && Length[rhs[[1]]] == 3 &&
            lhs[[1, 1]] === rhs[[1, 3]] && lhs[[1, 3]] === rhs[[1, 1]] &&
            lhs[[1, 2]] === rhs[[1, 2]],
        With[{da = Lookup[distToExit, lhs[[1, 1]], Infinity],
              dc = Lookup[distToExit, lhs[[1, 3]], Infinity]},
            Which[da > dc, rhs, da < dc, lhs, True, $Failed]
        ],
        True, $Failed
    ]

ResolveOneOr[_, _] := $Failed

ResolveOrByGraphDistance[Eqs_Association, triple:{_, _, _}] :=
    Module[{graph, exitVerts, distToExit, orTerms, resolved = {},
            unresolved = {}, newEE, newOR},
        If[triple[[3]] === True, Return[triple, Module]];
        graph = Lookup[Eqs, "auxiliaryGraph", None];
        exitVerts = Lookup[Eqs, "auxExitVertices", {}];
        If[graph === None || exitVerts === {},
            MFGPrint["ResolveOrByGraphDistance: no graph data available"];
            Return[triple, Module]
        ];
        distToExit = Association @ Table[
            v -> Min[GraphDistance[graph, v, #]& /@ exitVerts],
            {v, VertexList[graph]}
        ];
        orTerms = If[Head[triple[[3]]] === And,
            List @@ triple[[3]], {triple[[3]]}];
        Do[
            With[{res = ResolveOneOr[term, distToExit]},
                If[res =!= $Failed,
                    AppendTo[resolved, res],
                    AppendTo[unresolved, term]
                ]
            ],
            {term, orTerms}
        ];
        MFGPrint["ResolveOrByGraphDistance: resolved ", Length[resolved],
            " of ", Length[orTerms], " Or-conditions"];
        newEE = If[resolved === {},
            triple[[1]],
            And @@ Join[
                If[triple[[1]] === True, {},
                    If[Head[triple[[1]]] === And, List @@ triple[[1]],
                        {triple[[1]]}]],
                resolved
            ]
        ];
        newOR = If[unresolved === {}, True, And @@ unresolved];
        {newEE, triple[[2]], newOR}
    ]

(* --- DataToEquations: main converter --- *)

DataToEquations[Data_Association] :=
    Module[{verticesList, adjacencyMatrix, entryVerticesFlows, exitVerticesCosts,
         switchingCosts, graph, entryVertices, auxEntryVertices, exitVertices,
         auxExitVertices, entryEdges, exitEdges, auxiliaryGraph, auxEdgeList,
         edgeList, auxVerticesList, auxPairs, auxTriples, EntryDataAssociation,
         ExitCosts, js, us, jts, SignedFlows, SwitchingCosts, IneqJs, IneqJts,
         AltFlows, AltTransitionFlows, splittingPairs, EqBalanceSplittingFlows,
         BalanceSplittingFlows, NoDeadStarts, RuleBalanceGatheringFlows, BalanceGatheringFlows,
         EqBalanceGatheringFlows, EqEntryIn, RuleEntryValues, RuleEntryOut, RuleExitFlowsIn,
         RuleExitValues, EqValueAuxiliaryEdges, IneqSwitchingByVertex, AltOptCond,
         Nlhs, ModuleVars, ModuleVarsNames, Nrhs, costpluscurrents, EqGeneral,
         inAuxEntryPairs, outAuxEntryPairs, inAuxExitPairs, outAuxExitPairs, 
        pairs, halfPairs, consistentCosts},
        verticesList = Lookup[Data, "Vertices List", {}];
        adjacencyMatrix = Lookup[Data, "Adjacency Matrix", {}];
        entryVerticesFlows = Lookup[Data, "Entrance Vertices and Flows",
             {}];
        exitVerticesCosts = Lookup[Data, "Exit Vertices and Terminal Costs",
             {}];
        switchingCosts = Lookup[Data, "Switching Costs", {}];
        (* Graph construction *)
        graph = AdjacencyGraph[verticesList, adjacencyMatrix, VertexLabels
             -> "Name", DirectedEdges -> False];
        entryVertices = First /@ entryVerticesFlows;
        exitVertices = First /@ exitVerticesCosts;
(* Define auxiliary vertices for the entry and exit vertices 
    *)
        auxEntryVertices = Symbol["en" <> ToString[#]]& /@ entryVertices
            ;
        auxExitVertices = Symbol["ex" <> ToString[#]]& /@ exitVertices
            ;
        entryEdges = MapThread[UndirectedEdge, {auxEntryVertices, entryVertices
            }];
        exitEdges = MapThread[UndirectedEdge, {exitVertices, auxExitVertices
            }];
        auxiliaryGraph = EdgeAdd[graph, Join[entryEdges, exitEdges]];
            
        auxEdgeList = EdgeList[auxiliaryGraph];
        edgeList = EdgeList[graph];
        auxVerticesList = VertexList[auxiliaryGraph];
(* Arguments for the relevant functions: flows, values, and transition flows 
    
    
    
    *)
        halfPairs = List @@@ edgeList;
        inAuxEntryPairs = List @@@ entryEdges;
        outAuxExitPairs = List @@@ exitEdges;
        inAuxExitPairs = Reverse /@ outAuxExitPairs;
        outAuxEntryPairs = Reverse /@ inAuxEntryPairs;
        pairs = Join[halfPairs, Reverse /@ halfPairs];
        auxPairs = Join[inAuxEntryPairs, outAuxEntryPairs, inAuxExitPairs,
             outAuxExitPairs, pairs];
        (* Insert the vertex in pairs of adjacent vertices *)
        auxTriples = Flatten[Insert[#, 2] /@ Permutations[AdjacencyList[
            auxiliaryGraph, #], {2}]& /@ auxVerticesList, 1];
        (* Prepare boundary data *)
        EntryDataAssociation = RoundValues @ AssociationThread[inAuxEntryPairs,
             Last /@ entryVerticesFlows];
        ExitCosts = AssociationThread[auxExitVertices, Last /@ exitVerticesCosts
            ];
        (* Variables *)
        js = j[Sequence @@ #]& /@ auxPairs;
        jts = j[Sequence @@ #]& /@ auxTriples;
        us = u[Sequence @@ #]& /@ auxPairs;
        (* Signed flows for "half" of the variables *)
        SignedFlows = AssociationMap[j @@ # - j @@ Reverse @ #&, Join[
            inAuxEntryPairs, outAuxExitPairs, halfPairs]];
        SwitchingCosts = AssociationMap[0&, auxTriples];
        If[switchingCosts =!= {},
            AssociateTo[SwitchingCosts, AssociationThread[Most /@ switchingCosts,
                 Last /@ switchingCosts]]
        ];
        consistentCosts = IsSwitchingCostConsistent[Normal @ SwitchingCosts
            ];
        Which[
            consistentCosts === False,
                Message[DataToEquations::switchingcosts];
                Return[$Failed, Module]
            ,
            consistentCosts =!= True,
                MFGPrint["Switching costs conditions are ", consistentCosts
                    ]
        ];
        IneqJs = And @@ (# >= 0& /@ js);
        IneqJts = And @@ (# >= 0& /@ jts);
        (* Use only one representative per symmetric pair/triple to avoid
           generating duplicate conditions: AltFlowOp[j][{a,b}] and
           AltFlowOp[j][{b,a}] produce identical Or-conditions. *)
        AltFlows = And @@ (AltFlowOp[j] /@ Join[
            inAuxEntryPairs, outAuxExitPairs, halfPairs]);
        AltTransitionFlows = And @@ (AltFlowOp[j] /@
            Select[auxTriples, OrderedQ[{First[#], Last[#]}]&]);
        splittingPairs = Join[inAuxEntryPairs, inAuxExitPairs, pairs]
            ;
        BalanceSplittingFlows = (j @@ # - Total[j @@@ FlowSplitting[auxTriples
            ][#]])& /@ splittingPairs;
        EqBalanceSplittingFlows = Simplify /@ (And @@ ((# == 0)& /@ BalanceSplittingFlows
            ));
        NoDeadStarts = Join[outAuxEntryPairs, outAuxExitPairs, pairs]
            ;
        RuleBalanceGatheringFlows = Association[(j @@ # -> Total[j @@@
             FlowGathering[auxTriples][#]])& /@ NoDeadStarts];
        (* Equations for the exit currents at the entry vertices *)
        BalanceGatheringFlows = ((-j @@ # + Total[j @@@ FlowGathering[
            auxTriples][#]])& /@ NoDeadStarts);
        EqBalanceGatheringFlows = Simplify /@ (And @@ (# == 0& /@ BalanceGatheringFlows
            ));
        (* Incoming currents *)
        EqEntryIn = (j @@ # == EntryDataAssociation[#])& /@ inAuxEntryPairs
            ;
(* Outgoing flows at entrances and incoming flows at exits are zero 
    *)
        RuleEntryOut = Association[(j @@ # -> 0)& /@ outAuxEntryPairs
            ];
        RuleExitFlowsIn = Association[(j @@ # -> 0)& /@ inAuxExitPairs
            ];
        (* Exit values at exit vertices *)
        RuleExitValues = AssociationThread[u @@@ (Reverse /@ outAuxExitPairs
            ), Last /@ exitVerticesCosts];
        RuleExitValues = Join[RuleExitValues, AssociationThread[u @@@
             outAuxExitPairs, Last /@ exitVerticesCosts]];
        RuleEntryValues = AssociationThread[u @@@ outAuxEntryPairs, u
             @@@ inAuxEntryPairs];
        EqValueAuxiliaryEdges = And @@ ((u @@ # == u @@ Reverse[#])& 
            /@ Join[inAuxEntryPairs, outAuxExitPairs]);
        IneqSwitchingByVertex = IneqSwitch[u, SwitchingCosts] @@@ TransitionsAt[
            auxiliaryGraph, #]& /@ verticesList;
        IneqSwitchingByVertex = And @@@ IneqSwitchingByVertex;
        AltOptCond = And @@ AltSwitch[j, u, SwitchingCosts] @@@ auxTriples
            ;
        Nlhs = Flatten[u @@ # - u @@ Reverse @ # + SignedFlows[#]& /@
             halfPairs];
        Nrhs = Flatten[SignedFlows[#] - Sign[SignedFlows[#]] Cost[SignedFlows[
            #], #]& /@ halfPairs];
        (* Cost-plus-currents for the general (non-critical) case *)
        costpluscurrents = Table[Symbol["cpc" <> ToString[k]], {k, 1,
             Length @ edgeList}];
        EqGeneral = And @@ (MapThread[Equal, {Nlhs, costpluscurrents}
            ]);
        costpluscurrents = AssociationThread[costpluscurrents, Nrhs];
            
        (* Build output association *)
        ModuleVars = {graph, pairs, entryVertices, auxEntryVertices, 
            exitVertices, auxExitVertices, entryEdges, exitEdges, auxiliaryGraph,
             auxEdgeList, edgeList, auxVerticesList, auxTriples, EntryDataAssociation,
             ExitCosts, js, us, jts, SignedFlows, SwitchingCosts, IneqJs, IneqJts,
             AltFlows, AltTransitionFlows, splittingPairs, EqBalanceSplittingFlows,
             BalanceSplittingFlows, NoDeadStarts, RuleBalanceGatheringFlows, BalanceGatheringFlows,
             EqBalanceGatheringFlows, EqEntryIn, RuleEntryOut, RuleEntryValues, RuleExitFlowsIn,
             RuleExitValues, EqValueAuxiliaryEdges, IneqSwitchingByVertex, AltOptCond,
             Nlhs, Nrhs, costpluscurrents, EqGeneral};
        ModuleVarsNames = {"graph", "pairs", "entryVertices", "auxEntryVertices",
             "exitVertices", "auxExitVertices", "entryEdges", "exitEdges", "auxiliaryGraph",
             "auxEdgeList", "edgeList", "auxVerticesList", "auxTriples", "EntryDataAssociation",
             "ExitCosts", "js", "us", "jts", "SignedFlows", "SwitchingCosts", "IneqJs",
             "IneqJts", "AltFlows", "AltTransitionFlows", "splittingPairs", "EqBalanceSplittingFlows",
             "BalanceSplittingFlows", "NoDeadStarts", "RuleBalanceGatheringFlows",
             "BalanceGatheringFlows", "EqBalanceGatheringFlows", "EqEntryIn", "RuleEntryOut",
             "RuleEntryValues", "RuleExitFlowsIn", "RuleExitValues", "EqValueAuxiliaryEdges",
             "IneqSwitchingByVertex", "AltOptCond", "Nlhs", "Nrhs", "costpluscurrents",
             "EqGeneral"};
        Join[Data, AssociationThread[ModuleVarsNames, ModuleVars]]
    ];

(* --- Utility --- *)

NumberVectorQ[j_] :=
    VectorQ[j, NumberQ];

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
    Module[{Kirchhoff, EqEntryIn = Lookup[Eqs, "EqEntryIn", True], BalanceGatheringFlows
         = Lookup[Eqs, "BalanceGatheringFlows", {}], BalanceSplittingFlows = 
        Lookup[Eqs, "BalanceSplittingFlows", {}], RuleExitFlowsIn = Lookup[Eqs,
         "RuleExitFlowsIn", {}], RuleEntryOut = Lookup[Eqs, "RuleEntryOut", {
        }], BM, KM, vars, cost, CCost},
        (* 
           Construct the linear system ensuring flow conservation.
           Clever Cancellation: By summing BalanceGatheringFlows and BalanceSplittingFlows 
           for each edge, the internal edge-flow variables j(u,v) cancel out algebraically.
           This leaves a system of equations purely in terms of transition flows j(u,v,w).
           
           This reduction is critical for the Hessian Riemannian Flow (HRF) solver, 
           as it provides a minimal, independent coordinate basis (the state vector ϑ 
           in the HRF paper) while ensuring Kirchhoff's law is satisfied by construction.
        *)
        Kirchhoff = Join[EqEntryIn, (# == 0& /@ (BalanceGatheringFlows
             + BalanceSplittingFlows))];
        Kirchhoff = Kirchhoff /. Join[RuleExitFlowsIn, RuleEntryOut];
            
        (* Extract all flow variables from the Kirchhoff system *)
        vars = Select[Variables[Kirchhoff /. Equal -> List], MatchQ[#,
             j[_, _, _] | j[_, _]]&];
        If[Length[vars] === 0,
            (* Degenerate case: no flow variables *)
            {0, 0, {}}
            ,
            {BM, KM} = CoefficientArrays[Kirchhoff, vars];
            {-BM, KM, vars}
        ]
    ];

(*cost, below, is a function that should be defined by the problem*)

GetKirchhoffMatrix[Eqs_] :=
    Module[{B, KM, vars, cost, CCost},
        {B, KM, vars} = GetKirchhoffLinearSystem[Eqs];
        If[Length[vars] === 0,
            Return[{B, KM, <||>, vars}, Module]
        ];
(* The third slot is a legacy zero-cost placeholder retained only for compatibility. 
    
    
    
    *)
        cost = AssociationThread[vars, (0&) /@ vars];
        CCost = cost /@ vars /. MapThread[Rule, {vars, #}]&;
        {B, KM, CCost, vars}
    ];

(* --- MFGPreprocessing --- *)

MFGPreprocessing[Eqs_] :=
    Module[{InitRules, RuleBalanceGatheringFlows, EqEntryIn, RuleEntryOut,
         RuleExitFlowsIn, RuleExitValues, EqValueAuxiliaryEdges, IneqSwitchingByVertex,
         AltOptCond, EqBalanceSplittingFlows, AltFlows, AltTransitionFlows, IneqJs,
         IneqJts, ModuleVarsNames, ModulesVars, NewSystem, Rules, EqGeneral,
        RuleEntryValues, temp, time},
        RuleBalanceGatheringFlows = Lookup[Eqs, "RuleBalanceGatheringFlows",
             $Failed];
        EqEntryIn = Lookup[Eqs, "EqEntryIn", $Failed];
        RuleEntryOut = Lookup[Eqs, "RuleEntryOut", $Failed];
        RuleExitFlowsIn = Lookup[Eqs, "RuleExitFlowsIn", $Failed];
        RuleExitValues = Lookup[Eqs, "RuleExitValues", $Failed];
        RuleEntryValues = Lookup[Eqs, "RuleEntryValues", $Failed];
        AltOptCond = Lookup[Eqs, "AltOptCond", $Failed];
        AltFlows = Lookup[Eqs, "AltFlows", $Failed];
        AltTransitionFlows = Lookup[Eqs, "AltTransitionFlows", $Failed
            ];
        IneqJs = Lookup[Eqs, "IneqJs", $Failed];
        IneqJts = Lookup[Eqs, "IneqJts", $Failed];
        EqGeneral = Lookup[Eqs, "EqGeneral", $Failed];
        IneqSwitchingByVertex = Lookup[Eqs, "IneqSwitchingByVertex", 
            $Failed];
        EqBalanceSplittingFlows = Lookup[Eqs, "EqBalanceSplittingFlows",
             $Failed];
        EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", 
            $Failed];
        MFGPrint[Eqs["auxiliaryGraph"]];
        (* First rules: entry currents *)
        InitRules = Association[Flatten[ToRules /@ EqEntryIn]];
        AssociateTo[InitRules, RuleEntryOut];
        AssociateTo[InitRules, RuleExitFlowsIn];
        AssociateTo[InitRules, RuleEntryValues];
        AssociateTo[InitRules, RuleBalanceGatheringFlows];
        AssociateTo[InitRules, RuleExitValues];
        (* When all switching costs are zero, IneqSwitchingByVertex
           reduces to pairwise u[a,b]<=u[c,b] && u[c,b]<=u[a,b],
           i.e. equalities. Skip expensive Simplify and extract directly. *)
        If[AllTrue[Values[Lookup[Eqs, "SwitchingCosts", <||>]], # === 0 &],
            {time, IneqSwitchingByVertex} = AbsoluteTiming[
                Module[{items = IneqSwitchingByVertex /. InitRules, ineqs, eqs},
                    ineqs = Flatten[If[Head[#] === And, List @@ #, {#}]& /@ items];
                    eqs = DeleteDuplicatesBy[
                        Cases[ineqs, LessEqual[a_, b_] :> (a == b)],
                        Sort @* (List @@ #&)
                    ];
                    If[eqs === {}, True, And @@ eqs]
                ]
            ];
            MFGPrint["Switching costs (zero): extracted equalities in ",
                time, " seconds."];
            ,
            temp = MFGPrintTemporary["Simplifying inequalities involving Switching Costs...",
                 IneqSwitchingByVertex];
            With[{items = IneqSwitchingByVertex /. InitRules},
                {time, IneqSwitchingByVertex} = AbsoluteTiming[And @@ MFGParallelMap[
                    Simplify, items]]
            ];
            NotebookDelete[temp];
            MFGPrint["Switching costs simplified in ", time, " seconds."];
        ];
        EqBalanceSplittingFlows = EqBalanceSplittingFlows /. InitRules
            ;
        EqValueAuxiliaryEdges = EqValueAuxiliaryEdges /. InitRules;
        temp = MFGPrintTemporary["Solving some balance equations: ", 
            EqBalanceSplittingFlows && EqValueAuxiliaryEdges];
        {time, Rules} = AbsoluteTiming[First @ Solve[EqBalanceSplittingFlows
             && EqValueAuxiliaryEdges] // Quiet];
        InitRules = Join[InitRules /. Rules, Association @ Rules];
        NotebookDelete[temp];
        MFGPrint["Balance equations solved in ", time, " seconds."];
        With[{ineqTriple = SystemToTriple[IneqSwitchingByVertex]},
            NewSystem = {ineqTriple[[1]] && EqGeneral,
                IneqJts && IneqJs && ineqTriple[[2]],
                AltFlows && AltTransitionFlows && AltOptCond && ineqTriple[[3]]}
        ];
        (* Resolve Or-conditions using graph distance before TripleClean
           substitutes variables and obscures the j[a,b]==0 patterns *)
        {time, NewSystem} = AbsoluteTiming[
            ResolveOrByGraphDistance[Eqs, NewSystem]];
        MFGPrint["Graph-distance Or-resolution took ", time, " seconds."];
        temp = MFGPrintTemporary["TripleClean will work on\n", NewSystem,
             "\nwith: ", InitRules];
        {time, {NewSystem, InitRules}} = AbsoluteTiming @ TripleClean[
            {NewSystem, InitRules}];
        NotebookDelete[temp];
        MFGPrint["TripleClean took ", time, " seconds."];
        ModuleVarsNames = {"InitRules", "NewSystem"};
        ModulesVars = {InitRules, NewSystem};
        Join[Eqs, AssociationThread[ModuleVarsNames, ModulesVars]]
    ];

(* --- Direct numerical solver for zero-switching-cost critical congestion --- *)
(* When all switching costs are zero and Ncpc=0 (critical congestion at zero flow):
   - All values at a vertex are equal: u[a,b] = V(b)
   - EqGeneral: V(b) - V(a) = j[a,b] - j[b,a]
   - Complementarity: j[a,b]*j[b,a] = 0
   - Flow direction determined by graph distance to exit
   This reduces to a linear system solvable in O(V^2). *)

DirectCriticalSolver[Eqs_Association] :=
    Module[{graph, exitVerts, entryVerts, auxPairs, auxTriples,
            us, js, jts, distToExit,
            RuleEntryOut, RuleExitFlowsIn, RuleEntryValues, RuleExitValues,
            RuleBalanceGatheringFlows, EqGeneral, EqEntryIn,
            EqBalanceSplittingFlows, EqValueAuxiliaryEdges,
            AltFlows, AltTransitionFlows,
            zeroRules, eqSystem, allVars, sol, result},
        graph = Lookup[Eqs, "auxiliaryGraph", None];
        exitVerts = Lookup[Eqs, "auxExitVertices", {}];
        entryVerts = Lookup[Eqs, "auxEntryVertices", {}];
        us = Lookup[Eqs, "us", {}];
        js = Lookup[Eqs, "js", {}];
        jts = Lookup[Eqs, "jts", {}];
        auxPairs = Cases[us, u[a_, b_] :> {a, b}];
        auxTriples = Cases[jts, j[a_, b_, c_] :> {a, b, c}];
        RuleEntryOut = Lookup[Eqs, "RuleEntryOut", {}];
        RuleExitFlowsIn = Lookup[Eqs, "RuleExitFlowsIn", {}];
        RuleEntryValues = Lookup[Eqs, "RuleEntryValues", <||>];
        RuleExitValues = Lookup[Eqs, "RuleExitValues", <||>];
        RuleBalanceGatheringFlows = Lookup[Eqs, "RuleBalanceGatheringFlows", <||>];
        EqGeneral = Lookup[Eqs, "EqGeneral", True];
        EqEntryIn = Lookup[Eqs, "EqEntryIn", True];
        EqBalanceSplittingFlows = Lookup[Eqs, "EqBalanceSplittingFlows", True];
        EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", True];

        (* Phase 1: Compute vertex distances to exit *)
        distToExit = Association @ Table[
            v -> Min[GraphDistance[graph, v, #]& /@ exitVerts],
            {v, VertexList[graph]}
        ];

        (* Phase 2: Determine which flows are zero from graph distance.
           For edge {a,b}: if dist(a) < dist(b), then j[a,b]=0 (backward).
           For equidistant edges, both directions may carry flow.
           Transition flows are left free — determined by Kirchhoff. *)
        zeroRules = Association[];
        Do[With[{a = pair[[1]], b = pair[[2]]},
            With[{da = Lookup[distToExit, a, Infinity],
                  db = Lookup[distToExit, b, Infinity]},
                If[da < db, zeroRules[j[a, b]] = 0]
            ]], {pair, auxPairs}];
        MFGPrint["DirectCriticalSolver: ", Length[zeroRules],
            " flows set to zero by graph distance"];

        (* Phase 3: Build the full equation system.
           Substitute: Ncpc=0 (critical congestion at zero flow),
           known boundary rules, zero switching cost equalities
           (u[a,b]=u[c,b] at each vertex), and direction zeros. *)
        eqSystem = And @@ Flatten[{
            (* EqGeneral at Ncpc=0: j[a,b]-j[b,a]+u[a,b]-u[b,a]=0 *)
            List @@ (EqGeneral /. Thread[
                Cases[Keys[Lookup[Eqs, "costpluscurrents", <||>]], _] -> 0]),
            (* Entry flow *)
            If[Head[EqEntryIn] === And, List @@ EqEntryIn, {EqEntryIn}],
            (* Balance/splitting *)
            If[Head[EqBalanceSplittingFlows] === And,
                List @@ EqBalanceSplittingFlows,
                If[EqBalanceSplittingFlows === True, {}, {EqBalanceSplittingFlows}]],
            (* Value auxiliary edges *)
            If[Head[EqValueAuxiliaryEdges] === And,
                List @@ EqValueAuxiliaryEdges,
                If[EqValueAuxiliaryEdges === True, {}, {EqValueAuxiliaryEdges}]],
            (* Boundary values *)
            KeyValueMap[Equal, RuleExitValues],
            KeyValueMap[Equal, RuleEntryValues],
            (* Entry/exit flow rules *)
            (Equal @@@ Normal[RuleEntryOut]),
            (Equal @@@ Normal[RuleExitFlowsIn]),
            (* Gathering flow rules *)
            (Equal @@@ Normal[RuleBalanceGatheringFlows]),
            (* Direction zeros *)
            (Equal[#, 0]& /@ Keys[zeroRules]),
            (* Zero switching costs: u[a,b] = u[c,b] at each vertex *)
            (* Group by second element and equate all *)
            Module[{byVertex = GroupBy[auxPairs, Last]},
                Flatten @ KeyValueMap[
                    Function[{v, pairs},
                        If[Length[pairs] > 1,
                            (u @@ # == u @@ pairs[[1]])& /@ Rest[pairs],
                            {}
                        ]
                    ],
                    byVertex
                ]
            ]
        }];

        (* Phase 4: Solve with non-negativity on all flow variables *)
        allVars = Join[us, js, jts];
        With[{flowVars = Join[js, jts]},
            With[{ineqs = And @@ ((# >= 0 &) /@ flowVars)},
                With[{inst = Quiet @ FindInstance[eqSystem && ineqs, allVars, Reals]},
                    If[inst === {} || !MatchQ[inst, {{__Rule}, ___}],
                        MFGPrint["DirectCriticalSolver: FindInstance failed"];
                        Return[$Failed, Module]
                    ];
                    result = Association @ First @ inst;
                ]
            ]
        ];
        Do[If[!KeyExistsQ[result, var], result[var] = 0], {var, allVars}];
        (* Only resolve transitive chains if there are symbolic values *)
        If[!AllTrue[Values[result], NumericQ],
            result = Expand /@ FixedPoint[Function[r, ReplaceAll[r] /@ r], result, 10]
        ];
        result = Join[KeyTake[result, us], KeyTake[result, js], KeyTake[result, jts]];
        result
    ]

(* --- CriticalCongestionSolver --- *)

CriticalCongestionSolver[$Failed, ___] :=
    $Failed

CriticalCongestionSolver[Eqs_] :=
    Module[{PreEqs, js, AssoCritical, time, temp, status,
         resultKind, message, solution, comparisonData},
        ClearSolveCache[];
        (* Try direct solver for zero-switching-cost, fully numeric networks *)
        If[AllTrue[Values[Lookup[Eqs, "SwitchingCosts", <|_ -> 1|>]], # === 0 &] &&
            Lookup[Eqs, "auxiliaryGraph", None] =!= None &&
            AllTrue[Flatten[Lookup[Eqs, "Entrance Vertices and Flows", {}]], NumericQ] &&
            AllTrue[Flatten[Lookup[Eqs, "Exit Vertices and Terminal Costs", {}]], NumericQ],
            MFGPrint["Attempting direct critical solver (zero switching costs)..."];
            {time, AssoCritical} = AbsoluteTiming[DirectCriticalSolver[Eqs]];
            If[AssociationQ[AssoCritical] && AssoCritical =!= $Failed,
                MFGPrint["Direct critical solver completed in ", time, " seconds."];
                status = CheckFlowFeasibility[AssoCritical];
                If[status === "Feasible",
                    (* Run MFGPreprocessing so downstream solvers (NonLinearSolver)
                       have the symbolic InitRules, NewSystem, costpluscurrents that
                       MFGSystemSolver expects to specialize per-iteration. *)
                    If[KeyExistsQ[Eqs, "InitRules"],
                        PreEqs = Eqs,
                        {time, PreEqs} = AbsoluteTiming @ MFGPreprocessing[Eqs];
                        MFGPrint["Preprocessing (for downstream solvers) took ", time, " seconds."]
                    ];
                    solution = AssoCritical;
                    comparisonData = BuildSolverComparisonData[PreEqs, solution];
                    Return[
                        Join[PreEqs, MakeSolverResult["CriticalCongestion", "Success",
                            status, None, solution,
                            Join[comparisonData, <|"AssoCritical" -> AssoCritical,
                                "Status" -> status|>]]],
                        Module
                    ],
                    MFGPrint["Direct solver produced infeasible result, falling back to symbolic solver."]
                ],
                MFGPrint["Direct solver failed, falling back to symbolic solver."]
            ]
        ];
        (* Standard symbolic solver path *)
        If[KeyExistsQ[Eqs, "InitRules"],
            PreEqs = Eqs
            ,
            temp = MFGPrintTemporary["Preprocessing..."];
            {time, PreEqs} = AbsoluteTiming @ MFGPreprocessing[Eqs];
            NotebookDelete[temp];
            MFGPrint["Preprocessing took ", time, " seconds to terminate."
                ];
        ];
        js = Lookup[PreEqs, "js", $Failed];
        AssoCritical = MFGSystemSolver[PreEqs][AssociationThread[js,
            0 js]];
        status = CheckFlowFeasibility[AssoCritical];
        resultKind =
            If[AssoCritical === Null,
                "Failure"
                ,
                "Success"
            ];
        message =
            If[AssoCritical === Null,
                "NoSolution"
                ,
                None
            ];
        solution =
            If[AssociationQ[AssoCritical],
                AssoCritical
                ,
                Missing["NotAvailable"]
            ];
        comparisonData = BuildSolverComparisonData[PreEqs, solution];
        Join[PreEqs, MakeSolverResult["CriticalCongestion", resultKind,
             status, message, solution, Join[comparisonData, <|"AssoCritical" -> AssoCritical, "Status"
             -> status|>]]]
    ];

(* --- MFGSystemSolver --- *)

MFGSystemSolver[Eqs_][approxJs_] :=
    Module[{NewSystem, InitRules, pickOne, vars, System, Ncpc, costpluscurrents,
         us, js, jts, jjtsR, usR, time, temp, ineqsByTransition},
        ClearSolveCache[];
        us = Lookup[Eqs, "us", $Failed];
        js = Lookup[Eqs, "js", $Failed];
        jts = Lookup[Eqs, "jts", $Failed];
        InitRules = Lookup[Eqs, "InitRules", $Failed];
        NewSystem = Lookup[Eqs, "NewSystem", $Failed];
        costpluscurrents = Lookup[Eqs, "costpluscurrents", $Failed];
        {time, Ncpc} = AbsoluteTiming[RoundValues @ (Expand /@ (costpluscurrents
             /. approxJs))];
        MFGPrint["MFGSS: Calculated the cost plus currents for the flow in ",
             time, " seconds."];
        InitRules = Expand /@ (InitRules /. Ncpc);
        NewSystem = NewSystem /. Ncpc;
        (* Any False component makes the full {equalities, inequalities, alternatives}
           triple infeasible, so stop before attempting inequality bucketing. *)
        If[MemberQ[NewSystem, False],
            Message[MFGSystemSolver::nosolution];
            Return[Null, Module]
        ];
(* Retrieve some equalities from the inequalities: group by transition flow 
    
    
    
    *)
        temp = MFGPrintTemporary["MFGSS: Selecting inequalities by transition flow..."
            ];
        If[NewSystem[[2]] === True,
            ineqsByTransition = ConstantArray[True, Length[jts]];
            time = 0.
            ,
            {time, ineqsByTransition} =
                AbsoluteTiming[
                    Module[{ineqList, index},
                        ineqList = If[Head[NewSystem[[2]]] === And,
                            List @@ NewSystem[[2]], {NewSystem[[2]]}];
                        index = Association[# -> {} & /@ jts];
                        Do[
                            Do[
                                If[!FreeQ[ineq, jt],
                                    index[jt] = Append[index[jt], ineq]],
                                {jt, jts}
                            ],
                            {ineq, ineqList}
                        ];
                        (If[# === {}, True, And @@ #]&) /@ Values[index]
                    ]
                ]
        ];
        NotebookDelete[temp];
        MFGPrint["MFGSS: Selecting inequalities by transition flow took ",
             time, " seconds. ", Length[ineqsByTransition]];
        temp = MFGPrintTemporary["MFGSS: Simplifying inequalities by transition flow..."
            ];
        (* Pre-substitute known rules: many inequalities become True immediately *)
        ineqsByTransition = ineqsByTransition /. InitRules;
        ineqsByTransition = Replace[ineqsByTransition,
            expr_ /; expr =!= True :> Expand[expr], {1}];
        {time, ineqsByTransition} = AbsoluteTiming[DeduplicateByComplexity[
            MFGParallelMap[Simplify,
                Select[ineqsByTransition, # =!= True &]]]];
        (* Restore True entries for any that were already resolved *)
        NotebookDelete[temp];
        MFGPrint["MFGSS: Simplifying inequalities by transition flow took ",
             time, " seconds. "];
        NewSystem[[2]] = (And @@ ineqsByTransition);
        NewSystem = SystemToTriple[And @@ NewSystem];
        {NewSystem, InitRules} = TripleClean[{NewSystem, InitRules}];

        {NewSystem, InitRules} = DNFSolveStep[{NewSystem, InitRules}]
            ;
        MFGPrint["Simplifying system..."];
        System = Simplify[And @@ NewSystem];
        MFGPrint["Simplifying done. TripleClean..."];
        {System, InitRules} = TripleClean[{SystemToTriple[System], InitRules
            }];
        System = And @@ System;
        Which[
            System === False,
                Message[MFGSystemSolver::nosolution];
                Return[Null, Module]
            ,
            System =!= True,
                MFGPrint["MFGSS: (Possibly) Multiple solutions:\n", System
                    ];
                (* Extract all unsolved variables from System *)
                vars = Select[Variables[System], MatchQ[#, j[_, _, _]
                     | j[_, _] | u[_, _] | u[_, _, _]]&];
                vars = Complement[vars, Keys[InitRules]];
                jjtsR = Select[vars, MatchQ[#, j[_, _, _] | j[_, _]]&
                    ];
                usR = Select[vars, MatchQ[#, u[_, _, _] | u[_, _]]&];
                    
(* Pick one solution so that all the currents have numerical values 
    *)
                If[Length[vars] > 0,
                    (* Attempt to find a solution *)
                    pickOne = FindInstance[System && And @@ ((# > 0)&
                         /@ jjtsR), vars, Reals];
                    If[pickOne === {},
                        pickOne = FindInstance[System && And @@ ((# >=
                             0)& /@ jjtsR), vars, Reals]
                    ];
                    If[pickOne =!= {},
                        pickOne = Association @ First @ pickOne;
                        MFGPrint["MFGSS: Picked one solution: ", pickOne
                            ];
                        InitRules = Expand /@ Join[InitRules /. pickOne,
                             pickOne]
                        ,
                        MFGPrint["MFGSS: No feasible solution found"]
                            ;
                        Return[Null, Module]
                    ]
                    ,
(* All variables already solved — nothing to pick 
    *)
                    MFGPrint["MFGSS: All variables already determined in InitRules"
                        ]
                ]
            ,
            True,
                MFGPrint["MFGSS: System is ", System]
        ];
(* Resolve transitive rule chains, e.g. j[2,4] -> j[1,2,4]+j[3,2,4]
   with j[1,2,4] -> 50 becomes j[2,4] -> 50 *)
        InitRules =
            Expand /@
                FixedPoint[
                    Function[r,
                        ReplaceAll[r] /@ r
                    ]
                    ,
                    InitRules
                    ,
                    10
                ];
        InitRules = Join[KeyTake[InitRules, us], KeyTake[InitRules, js
            ], KeyTake[InitRules, jts]];
        InitRules
    ];

(* --- DNFSolveStep --- *)

DNFSolveStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, rules_}] :=
    {{EE, NN, OR}, rules}

DNFSolveStep[{{EE_, NN_, OO_}, rules_}] :=
    Module[{NewSystem, newrules, sorted = True, time, temp},
        {NewSystem, newrules} = TripleClean[{{EE, NN, OO}, rules}];
        If[NewSystem[[3]] === True,
            Return[{NewSystem, newrules}, Module]
            ,
            MFGPrint["Final: ", Length /@ NewSystem];
            sorted = DeduplicateByComplexity[NewSystem[[3]]];
            (* Sort Or-conditions by descending LeafCount: processing
               larger conditions first lets DNFReduce constrain the
               system earlier, improving early-exit rates. *)
            If[Head[sorted] === And,
                sorted = And @@ SortBy[List @@ sorted, -LeafCount[#] &]
            ];
        ];
        temp = MFGPrintTemporary["Final: Iterative DNF conversion on ",
             Length[sorted], " disjunctions..."];
        {time, NewSystem} = AbsoluteTiming[DNFReduce[And @@ Most[NewSystem
            ], sorted]];
        NotebookDelete[temp];
        MFGPrint["Final: Iterative DNF conversion on ", Length[sorted
            ], " disjunctions took ", time, " seconds to terminate."];
        {time, NewSystem} = AbsoluteTiming[ReduceDisjuncts[NewSystem]
            ];
        MFGPrint[
            "Final: ReduceDisjuncts returned "
            ,
            If[Head[NewSystem] === Or,
                Length[NewSystem]
                ,
                1
            ]
            ,
            " disjunct(s) in "
            ,
            time
            ,
            " seconds."
        ];
        NewSystem = BooleanConvert @ NewSystem;
        NewSystem = SystemToTriple[NewSystem];
        {NewSystem, newrules} = TripleClean[{NewSystem, newrules}];
        MFGPrint["Now: ", TimeObject[Now], " The new rules are: ", newrules,
             Length /@ NewSystem];
        {NewSystem, newrules}
    ];

(* --- SystemToTriple: decompose into equalities, inequalities, alternatives --- *)

SystemToTriple[True] = Table[True, 3]

SystemToTriple[False] = Table[False, 3]

SystemToTriple[system_] :=
    Which[
        Head[system] === And,
            Module[{args = List @@ system},
                {And @@ Cases[args, _Equal], And @@ DeleteCases[args,
                     _Equal | _Or], And @@ Cases[args, _Or]}
            ]
        ,
        Head[system] === Equal,
            {system, True, True}
        ,
        Head[system] === Or,
            {True, True, system}
        ,
        True,
            {True, system, True}
    ];

(* --- TripleStep / TripleClean: fixed-point simplification --- *)

(* SafeFirstSolve: returns a list of Rules from CachedSolve, or {} on failure *)

SafeFirstSolve[eq_] :=
    Module[{sol},
        sol = Quiet[CachedSolve[eq]];
        If[MatchQ[sol, {{__Rule}, ___}],
            First[sol]
            ,
            {}
        ]
    ];

TripleStep[{{EEs_, NNs_, ORs_}, rules_List}] :=
    TripleStep[{{EEs, NNs, ORs}, Association @ rules}]

TripleStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, rules_Association
    }] :=
    {{EE, NN, OR}, rules}

TripleStep[{{EE_, NN_?TrueQ, OR_?TrueQ}, rules_Association}] :=
    Module[{newrules = {}},
        newrules = SafeFirstSolve[EE /. rules];
        newrules = Join[rules /. newrules, Association @ newrules];
        {{True, NN, OR}, Expand /@ newrules}
    ];

TripleStep[{{True, NNs_, ORs_}, rules_Association}] :=
    {{True, NNs, ORs}, rules}

TripleStep[{{EEs_, NNs_, ORs_}, rules_Association}] :=
    Module[{EE = EEs /. rules, NN, OR, NNE, NNO, ORE, ORN, newrules =
         {}},
        If[EE =!= True && EE =!= False,
            newrules = SafeFirstSolve[EE]
        ];
        newrules = Expand /@ Join[rules /. newrules, Association @ newrules
            ];
        NN = Expand /@ (NNs /. newrules);
        OR = Expand /@ (ORs /. newrules);
        {NNE, NN, NNO} = SystemToTriple[NN];
        {ORE, ORN, OR} = SystemToTriple[OR];
        EE = NNE && ORE;
        {{EE, NN, OR}, newrules}
    ];

TripleClean[{{EE_, NN_, OR_}, rules_}] :=
    FixedPoint[TripleStep, {{EE, NN, OR}, rules}];

(* --- Backward compatibility aliases --- *)

Data2Equations = DataToEquations;

FinalStep = DNFSolveStep;

End[];
