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
payload key \"AssoCritical\". Options: \"ValidateSolution\" (default True), \
\"ValidationTolerance\" (default $CriticalSolverTolerance), and \"ValidationVerbose\" (default False).";

IsCriticalSolution::usage =
"IsCriticalSolution[Eqs] validates whether the critical-congestion solution stored \
in \"AssoCritical\" satisfies the full symbolic MFG constraint set (equalities, \
inequalities, alternatives/complementarity) and the critical EqGeneral residual. \
By default it returns True/False. Options: \"Tolerance\" (default $CriticalSolverTolerance), \
\"Verbose\" (default False), and \"ReturnReport\" (default False).";

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

(* --- Solver tolerance constant --- *)
(* Single source of truth for numeric tolerance across all critical-congestion
   solver paths: JFirst backend, coupled numeric backend, DirectCriticalSolver,
   and the IsCriticalSolution validation gate. The CriticalCongestionSolver
   "ValidationTolerance" option overrides this when set explicitly. *)
$CriticalSolverTolerance = 10^-6;

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

(* --- Numeric state compiler and adapters (internal) --- *)

ComputeDistanceToExitAssociation[graph_, exitVerts_List] :=
    Module[{verts, distMatrix},
        If[Head[graph] =!= Graph || exitVerts === {},
            Return[<||>, Module]
        ];
        verts = VertexList[graph];
        If[verts === {},
            Return[<||>, Module]
        ];
        distMatrix = Quiet @ Check[
            Normal @ GraphDistanceMatrix[graph, verts, exitVerts],
            $Failed
        ];
        If[distMatrix === $Failed || !MatrixQ[distMatrix],
            Return[<||>, Module]
        ];
        AssociationThread[verts, Min /@ distMatrix]
    ];

BuildSignedEdgeMatrixFromKirchhoff[Eqs_Association, kirchhoffVars_List] :=
    Module[{edgePairs, signedFlows, rules, signedExprs, placeholders, substituted, signedMatrix},
        edgePairs = List @@@ Lookup[Eqs, "edgeList", {}];
        If[edgePairs === {} || kirchhoffVars === {},
            Return[SparseArray[{}, {Length[edgePairs], Length[kirchhoffVars]}], Module]
        ];
        signedFlows = Lookup[Eqs, "SignedFlows", <||>];
        rules = Join[
            Lookup[Eqs, "RuleBalanceGatheringFlows", <||>],
            Lookup[Eqs, "RuleExitFlowsIn", <||>],
            Lookup[Eqs, "RuleEntryOut", <||>]
        ];
        signedExprs = (If[KeyExistsQ[signedFlows, #], signedFlows[#], 0] & /@ edgePairs) /. rules;
        placeholders = Table[Unique["sv"], {Length[kirchhoffVars]}];
        substituted = signedExprs /. Thread[kirchhoffVars -> placeholders];
        signedMatrix = Quiet @ Check[
            Last @ CoefficientArrays[substituted, placeholders],
            $Failed
        ];
        If[signedMatrix === $Failed,
            SparseArray[{}, {Length[edgePairs], Length[kirchhoffVars]}],
            If[Head[signedMatrix] === SparseArray,
                signedMatrix,
                SparseArray[signedMatrix]
            ]
        ]
    ];

ExtractDirectedEdgePairsFromAdjacencyMatrix[adjacency_, vertexList_List] :=
    Module[{n, rules, densePairs, edgePairs},
        If[vertexList === {},
            Return[{}, Module]
        ];
        n = Length[vertexList];
        edgePairs =
            Which[
                Head[adjacency] === SparseArray,
                    rules = Most[ArrayRules[adjacency]];
                    Cases[
                        rules,
                        ({i_Integer, j_Integer} -> val_) /;
                            1 <= i <= n && 1 <= j <= n && !PossibleZeroQ[val] :>
                            {vertexList[[i]], vertexList[[j]]}
                    ],
                MatrixQ[adjacency],
                    densePairs = Reap[
                        Do[
                            If[
                                !PossibleZeroQ[adjacency[[i, j]]],
                                Sow[{vertexList[[i]], vertexList[[j]]}]
                            ],
                            {i, 1, Min[Length[adjacency], n]},
                            {j, 1, Min[Length[adjacency[[i]]], n]}
                        ]
                    ][[2]];
                    If[densePairs === {}, {}, First[densePairs]],
                True,
                    {}
            ];
        DeleteDuplicates[edgePairs]
    ];

BuildCriticalDecouplingPartition[Eqs_Association, numericState_Association] :=
    Module[{jVars, uVars, kirchhoffVars, kirchhoffIndices, km, bRaw, bVec,
      adjacency, vertexList, edgePairs, directedGraph, cycleSample, isDAG,
      nodeOrder, nodePosition, edgeOrder, jVarToIndex, jVarOrder, cycleWitness},
        jVars = Lookup[numericState, "JVars", {}];
        uVars = Lookup[numericState, "UVars", {}];
        kirchhoffVars = Lookup[numericState, "KirchhoffVariables", {}];
        kirchhoffIndices = Lookup[numericState, "KirchhoffIndicesInFlowVector", {}];
        km = Lookup[numericState, "KirchhoffKM", SparseArray[{}, {0, 0}]];
        If[Head[km] =!= SparseArray,
            km = SparseArray[km]
        ];
        bRaw = Lookup[numericState, "KirchhoffB", {}];
        bVec = Developer`ToPackedArray @ N @ Which[
            Head[bRaw] === SparseArray, Normal[bRaw],
            VectorQ[bRaw], bRaw,
            True, {}
        ];

        adjacency = Lookup[Eqs, "Adjacency Matrix", {}];
        vertexList = Lookup[Eqs, "Vertices List", Lookup[numericState, "VertexList", {}]];
        edgePairs = ExtractDirectedEdgePairsFromAdjacencyMatrix[adjacency, vertexList];
        directedGraph = Graph[vertexList, DirectedEdge @@@ edgePairs, DirectedEdges -> True];
        cycleSample = Quiet @ Check[FindCycle[directedGraph, {1, Infinity}, 1], {}];
        isDAG = cycleSample === {};
        nodeOrder = If[isDAG, Quiet @ Check[TopologicalSort[directedGraph], {}], {}];
        If[!ListQ[nodeOrder], nodeOrder = {}];
        nodePosition = AssociationThread[nodeOrder, Range[Length[nodeOrder]]];
        edgeOrder =
            If[isDAG,
                SortBy[
                    edgePairs,
                    {
                        Lookup[nodePosition, First[#], Infinity] &,
                        Lookup[nodePosition, Last[#], Infinity] &,
                        ToString[First[#], InputForm] &,
                        ToString[Last[#], InputForm] &
                    }
                ],
                {}
            ];
        jVarToIndex = Lookup[numericState, "JVarToIndex",
            AssociationThread[jVars, Range[Length[jVars]]]
        ];
        jVarOrder =
            If[isDAG,
                Select[j @@@ edgeOrder, KeyExistsQ[jVarToIndex, #] &],
                {}
            ];
        cycleWitness =
            If[isDAG,
                Missing["NotApplicable"],
                If[ListQ[cycleSample] && cycleSample =!= {},
                    First[cycleSample],
                    Missing["NotAvailable"]
                ]
            ];

        <|
            "JVars" -> jVars,
            "UVars" -> uVars,
            "MassConstraints" -> <|
                "Aeq" -> km,
                "beq" -> bVec,
                "Variables" -> kirchhoffVars,
                "FlowVariableIndices" -> kirchhoffIndices,
                "LowerBounds" -> Developer`ToPackedArray @ ConstantArray[0., Length[kirchhoffVars]],
                "UpperBounds" -> Developer`ToPackedArray @ ConstantArray[Infinity, Length[kirchhoffVars]],
                "ConstraintKind" -> "KirchhoffOnly"
            |>,
            "TopologicalOrder" -> <|
                "IsDAG" -> isDAG,
                "NodeOrder" -> nodeOrder,
                "EdgeOrder" -> edgeOrder,
                "JVarOrder" -> jVarOrder,
                "CycleWitness" -> cycleWitness
            |>
        |>
    ];

BuildNumericState[Eqs_Association] :=
    Module[{graph, vertexList, edgeList, edgePairs, transitionTriples, js, jts, us,
      flowVariables, jVarToIndex, jtVarToIndex, uVarToIndex, flowVarToIndex,
      B, KM, kirchhoffVars, kirchhoffIndices, signedEdgeMatrix, statePairs,
      numericState},
        graph = Lookup[Eqs, "auxiliaryGraph", None];
        vertexList = If[Head[graph] === Graph, VertexList[graph], {}];
        edgeList = Lookup[Eqs, "edgeList", {}];
        edgePairs = List @@@ edgeList;
        transitionTriples = Lookup[Eqs, "auxTriples", {}];
        js = Lookup[Eqs, "js", {}];
        jts = Lookup[Eqs, "jts", {}];
        us = Lookup[Eqs, "us", {}];
        flowVariables = Join[js, jts];

        jVarToIndex = AssociationThread[js, Range[Length[js]]];
        jtVarToIndex = AssociationThread[jts, Range[Length[jts]]];
        uVarToIndex = AssociationThread[us, Range[Length[us]]];
        flowVarToIndex = AssociationThread[flowVariables, Range[Length[flowVariables]]];

        {B, KM, kirchhoffVars} = GetKirchhoffLinearSystem[Eqs];
        kirchhoffIndices = Lookup[flowVarToIndex, kirchhoffVars, Missing["NotAvailable"]];
        If[MemberQ[kirchhoffIndices, _Missing],
            kirchhoffIndices = {}
        ];

        signedEdgeMatrix = BuildSignedEdgeMatrixFromKirchhoff[Eqs, kirchhoffVars];
        statePairs = DeleteDuplicates @ Cases[us, u[a_, b_] :> {a, b}];

        numericState = <|
            "Version" -> 1,
            "VertexList" -> vertexList,
            "VertexIndex" -> AssociationThread[vertexList, Range[Length[vertexList]]],
            "EdgePairs" -> edgePairs,
            "EdgeIndex" -> AssociationThread[edgePairs, Range[Length[edgePairs]]],
            "TransitionTriples" -> transitionTriples,
            "TransitionIndex" -> AssociationThread[transitionTriples, Range[Length[transitionTriples]]],
            "StatePairs" -> statePairs,
            "StatePairIndex" -> AssociationThread[statePairs, Range[Length[statePairs]]],
            "JVars" -> js,
            "JVarToIndex" -> jVarToIndex,
            "JTVars" -> jts,
            "JTVarToIndex" -> jtVarToIndex,
            "UVars" -> us,
            "UVarToIndex" -> uVarToIndex,
            "FlowVariables" -> flowVariables,
            "FlowVarToIndex" -> flowVarToIndex,
            "KirchhoffVariables" -> kirchhoffVars,
            "KirchhoffVarToIndex" -> AssociationThread[kirchhoffVars, Range[Length[kirchhoffVars]]],
            "KirchhoffIndicesInFlowVector" -> kirchhoffIndices,
            "KirchhoffB" -> Developer`ToPackedArray @ N @ Which[
                Head[B] === SparseArray, Normal[B],
                VectorQ[B], B,
                True, {}
            ],
            "KirchhoffKM" ->
                If[Head[KM] === SparseArray,
                    KM,
                    SparseArray[KM]
                ],
            "SignedEdgeMatrix" -> signedEdgeMatrix,
            "DistanceToExit" -> ComputeDistanceToExitAssociation[graph, Lookup[Eqs, "auxExitVertices", {}]]
        |>;
        Join[
            numericState,
            <|"CriticalDecoupling" -> BuildCriticalDecouplingPartition[Eqs, numericState]|>
        ]
    ];

EncodeFlowAssociation[numericState_Association, assoc_Association] :=
    Module[{vars, vals},
        vars = Lookup[numericState, "FlowVariables", {}];
        If[vars === {},
            Return[Developer`ToPackedArray[{}], Module]
        ];
        vals = Lookup[assoc, vars, Missing["NotAvailable"]];
        If[VectorQ[vals, NumericQ],
            Developer`ToPackedArray @ N @ vals,
            Missing["NotAvailable"]
        ]
    ];

EncodeFlowAssociation[_, _] := Missing["NotAvailable"];

DecodeFlowVector[numericState_Association, flowVec_] :=
    Module[{vars, vals},
        vars = Lookup[numericState, "FlowVariables", {}];
        vals = Quiet @ Check[Developer`ToPackedArray @ N @ flowVec, $Failed];
        If[vals === $Failed || !VectorQ[vals, NumericQ] || Length[vals] =!= Length[vars],
            <||>,
            AssociationThread[vars, vals]
        ]
    ];

DecodeFlowVector[_, _] := <||>;

ComputeSignedEdgeFlowsFast[numericState_Association, flowVec_] :=
    Module[{indices, signedEdgeMatrix, kirchhoffVec},
        indices = Lookup[numericState, "KirchhoffIndicesInFlowVector", {}];
        signedEdgeMatrix = Lookup[numericState, "SignedEdgeMatrix", SparseArray[{}, {0, 0}]];
        If[Dimensions[signedEdgeMatrix][[1]] === 0,
            Return[Developer`ToPackedArray[{}], Module]
        ];
        If[indices === {} || !VectorQ[indices, IntegerQ],
            Return[Missing["NotAvailable"], Module]
        ];
        If[Length[flowVec] < Max[indices],
            Return[Missing["NotAvailable"], Module]
        ];
        kirchhoffVec = flowVec[[indices]];
        If[!VectorQ[kirchhoffVec, NumericQ],
            Return[Missing["NotAvailable"], Module]
        ];
        Developer`ToPackedArray @ N @ (signedEdgeMatrix . kirchhoffVec)
    ];

ComputeKirchhoffResidualFast[numericState_Association, flowVec_] :=
    Module[{indices, KM, B, kirchhoffVec, residual},
        indices = Lookup[numericState, "KirchhoffIndicesInFlowVector", {}];
        KM = Lookup[numericState, "KirchhoffKM", SparseArray[{}, {0, 0}]];
        B = Lookup[numericState, "KirchhoffB", Developer`ToPackedArray[{}]];
        If[indices === {} || !VectorQ[indices, IntegerQ],
            Return[Missing["NotAvailable"], Module]
        ];
        If[Length[flowVec] < Max[indices],
            Return[Missing["NotAvailable"], Module]
        ];
        kirchhoffVec = flowVec[[indices]];
        If[!VectorQ[kirchhoffVec, NumericQ],
            Return[Missing["NotAvailable"], Module]
        ];
        residual = Quiet @ Check[N @ Norm[KM . kirchhoffVec - B, Infinity], Missing["NotAvailable"]];
        If[NumericQ[residual], residual, Missing["NotAvailable"]]
    ];

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
        Module[{baseAssoc},
            baseAssoc = Join[Data, AssociationThread[ModuleVarsNames, ModuleVars]];
            Join[baseAssoc, <|"NumericState" -> BuildNumericState[baseAssoc]|>]
        ]
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

UtilityConjunctionTerms[expr_] :=
    Which[
        TrueQ[expr], {},
        Head[expr] === And, List @@ expr,
        True, {expr}
    ];

CanonicalUtilityOrderKey[expr_] := ToString[expr, InputForm];

BuildUtilityReductionData[Eqs_Association, initRules_Association, switchingExpr_, valueExpr_] :=
    Module[{switchingVals, uVars, rulesList, resolve, rulePairs, switchingPairs, valuePairs,
      utilityPairs, graph, classes, classData, representativeRules, representativeVariables,
      reducedEqGeneral, reducedSwitching, reducedValue},
        switchingVals = Values[Lookup[Eqs, "SwitchingCosts", <||>]];
        If[!AllTrue[switchingVals, # === 0 &],
            Return[
                <|
                    "Enabled" -> False,
                    "Reason" -> "NonZeroSwitchingCosts",
                    "UtilityClasses" -> {},
                    "RepresentativeRules" -> <||>,
                    "RepresentativeVariables" -> {},
                    "ClassCount" -> 0,
                    "ReducedVariableCount" -> 0,
                    "ReducedEqGeneral" -> Missing["NotAvailable"],
                    "ReducedSwitchingEqualities" -> Missing["NotAvailable"],
                    "ReducedEqValueAuxiliaryEdges" -> Missing["NotAvailable"]
                |>,
                Module
            ]
        ];
        uVars = Lookup[Eqs, "us", {}];
        If[uVars === {} || !AssociationQ[initRules],
            Return[
                <|
                    "Enabled" -> False,
                    "Reason" -> "MissingUtilityState",
                    "UtilityClasses" -> {},
                    "RepresentativeRules" -> <||>,
                    "RepresentativeVariables" -> {},
                    "ClassCount" -> 0,
                    "ReducedVariableCount" -> 0,
                    "ReducedEqGeneral" -> Missing["NotAvailable"],
                    "ReducedSwitchingEqualities" -> Missing["NotAvailable"],
                    "ReducedEqValueAuxiliaryEdges" -> Missing["NotAvailable"]
                |>,
                Module
            ]
        ];

        rulesList = Normal[initRules];
        resolve[expr_] := FixedPoint[
            Function[val, Quiet @ Check[val /. rulesList, val]],
            expr,
            12
        ];

        rulePairs = DeleteDuplicates @ Cases[
            rulesList,
            Rule[lhs_, rhs_] /;
                MatchQ[resolve[lhs], _u] &&
                MatchQ[resolve[rhs], _u] &&
                resolve[lhs] =!= resolve[rhs] :>
                    SortBy[{resolve[lhs], resolve[rhs]}, CanonicalUtilityOrderKey]
        ];

        switchingPairs = DeleteDuplicates @ Cases[
            UtilityConjunctionTerms[switchingExpr],
            Equal[lhs_, rhs_] /;
                MatchQ[resolve[lhs], _u] &&
                MatchQ[resolve[rhs], _u] &&
                resolve[lhs] =!= resolve[rhs] :>
                    SortBy[{resolve[lhs], resolve[rhs]}, CanonicalUtilityOrderKey]
        ];

        valuePairs = DeleteDuplicates @ Cases[
            UtilityConjunctionTerms[valueExpr],
            Equal[lhs_, rhs_] /;
                MatchQ[resolve[lhs], _u] &&
                MatchQ[resolve[rhs], _u] &&
                resolve[lhs] =!= resolve[rhs] :>
                    SortBy[{resolve[lhs], resolve[rhs]}, CanonicalUtilityOrderKey]
        ];

        utilityPairs = DeleteDuplicates @ Join[rulePairs, switchingPairs, valuePairs];
        graph = Graph[uVars, UndirectedEdge @@@ utilityPairs];
        classes = ConnectedComponents[graph];

        classData = Map[
            Function[class,
                Module[{resolvedClass, anchorValues, representative, classRules},
                    resolvedClass = DeleteDuplicates[resolve /@ class];
                    anchorValues = DeleteDuplicates @ Select[resolvedClass, FreeQ[#, _u] &];
                    representative =
                        If[anchorValues =!= {},
                            First @ SortBy[anchorValues, CanonicalUtilityOrderKey],
                            First @ SortBy[class, CanonicalUtilityOrderKey]
                        ];
                    classRules = Cases[class, var_ /; var =!= representative :> (var -> representative)];
                    <|
                        "Class" -> class,
                        "Representative" -> representative,
                        "Anchored" -> (anchorValues =!= {}),
                        "AnchorValues" -> anchorValues,
                        "Rules" -> classRules
                    |>
                ]
            ],
            classes
        ];

        representativeRules = Association @ Flatten[Lookup[classData, "Rules", {}]];
        representativeVariables = DeleteDuplicates @ Cases[
            Lookup[classData, "Representative", {}],
            _u
        ];
        reducedEqGeneral = Lookup[Eqs, "EqGeneral", True] /. Normal[representativeRules];
        reducedSwitching = switchingExpr /. Normal[representativeRules];
        reducedValue = valueExpr /. Normal[representativeRules];

        <|
            "Enabled" -> True,
            "Reason" -> None,
            "UtilityClasses" -> classes,
            "ClassData" -> classData,
            "RepresentativeRules" -> representativeRules,
            "RepresentativeVariables" -> representativeVariables,
            "ClassCount" -> Length[classes],
            "ReducedVariableCount" -> Length[representativeVariables],
            "ReducedEqGeneral" -> reducedEqGeneral,
            "ReducedSwitchingEqualities" -> reducedSwitching,
            "ReducedEqValueAuxiliaryEdges" -> reducedValue
        |>
    ];

(* --- MFGPreprocessing --- *)

MFGPreprocessing[Eqs_] :=
    Module[{InitRules, RuleBalanceGatheringFlows, EqEntryIn, RuleEntryOut,
         RuleExitFlowsIn, RuleExitValues, EqValueAuxiliaryEdges, IneqSwitchingByVertex,
         AltOptCond, EqBalanceSplittingFlows, AltFlows, AltTransitionFlows, IneqJs,
         IneqJts, ModuleVarsNames, ModulesVars, NewSystem, Rules, EqGeneral,
        utilityReductionData,
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
        utilityReductionData = BuildUtilityReductionData[
            Eqs,
            InitRules,
            IneqSwitchingByVertex,
            EqValueAuxiliaryEdges
        ];
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
        ModuleVarsNames = {"InitRules", "NewSystem", "UtilityReduction"};
        ModulesVars = {InitRules, NewSystem, utilityReductionData};
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

(* --- Critical numeric backend (internal, guarded) --- *)

$CriticalNumericBackendEnabled = False;

CriticalNumericBackendRequestedQ[Eqs_Association] :=
    Module[{mode},
        mode = Lookup[Eqs, "CriticalNumericBackendMode", Automatic];
        Which[
            mode === True || mode === "Force", True,
            mode === False || mode === "Disabled", False,
            mode === Automatic || mode === "Auto", TrueQ[$CriticalNumericBackendEnabled],
            True, False
        ]
    ];

CriticalNumericBackendEligibleQ[Eqs_Association] :=
    Module[{numericState, switchingVals, entryVals, exitVals, km, b, vars},
        numericState = Lookup[Eqs, "NumericState", Missing["NotAvailable"]];
        switchingVals = Values[Lookup[Eqs, "SwitchingCosts", <||>]];
        entryVals = Flatten[Lookup[Eqs, "Entrance Vertices and Flows", {}]];
        exitVals = Flatten[Lookup[Eqs, "Exit Vertices and Terminal Costs", {}]];
        km = Lookup[numericState, "KirchhoffKM", Missing["NotAvailable"]];
        b = Lookup[numericState, "KirchhoffB", Missing["NotAvailable"]];
        vars = Lookup[numericState, "KirchhoffVariables", {}];
        AssociationQ[numericState] &&
        AllTrue[switchingVals, # === 0 &] &&
        Lookup[Eqs, "auxiliaryGraph", None] =!= None &&
        AllTrue[entryVals, NumericQ] &&
        AllTrue[exitVals, NumericQ] &&
        Head[km] === SparseArray &&
        VectorQ[b, NumericQ] &&
        ListQ[vars] &&
        vars =!= {}
    ];

BuildCriticalLinearConstraints[Eqs_Association] :=
    Module[{graph, exitVerts, us, js, jts, auxPairs, distToExit, zeroRules, pair,
        RuleEntryOut, RuleExitFlowsIn, RuleEntryValues, RuleExitValues,
        RuleBalanceGatheringFlows, EqGeneral, EqEntryIn,
        EqBalanceSplittingFlows, EqValueAuxiliaryEdges,
        eqList, flowVars, allVars, ineqList},
        graph = Lookup[Eqs, "auxiliaryGraph", None];
        exitVerts = Lookup[Eqs, "auxExitVertices", {}];
        us = Lookup[Eqs, "us", {}];
        js = Lookup[Eqs, "js", {}];
        jts = Lookup[Eqs, "jts", {}];
        auxPairs = Cases[us, u[a_, b_] :> {a, b}];
        RuleEntryOut = Lookup[Eqs, "RuleEntryOut", {}];
        RuleExitFlowsIn = Lookup[Eqs, "RuleExitFlowsIn", {}];
        RuleEntryValues = Lookup[Eqs, "RuleEntryValues", <||>];
        RuleExitValues = Lookup[Eqs, "RuleExitValues", <||>];
        RuleBalanceGatheringFlows = Lookup[Eqs, "RuleBalanceGatheringFlows", <||>];
        EqGeneral = Lookup[Eqs, "EqGeneral", True];
        EqEntryIn = Lookup[Eqs, "EqEntryIn", True];
        EqBalanceSplittingFlows = Lookup[Eqs, "EqBalanceSplittingFlows", True];
        EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", True];

        distToExit = Lookup[
            Lookup[Eqs, "NumericState", <||>],
            "DistanceToExit",
            <||>
        ];
        If[!AssociationQ[distToExit] || distToExit === <||>,
            If[Head[graph] =!= Graph || exitVerts === {},
                Return[$Failed, Module]
            ];
            distToExit = Association @ Table[
                v -> Min[GraphDistance[graph, v, #] & /@ exitVerts],
                {v, VertexList[graph]}
            ];
        ];

        zeroRules = Association[];
        Do[
            pair = auxPairs[[k]];
            If[
                Lookup[distToExit, pair[[1]], Infinity] <
                    Lookup[distToExit, pair[[2]], Infinity],
                zeroRules[j @@ pair] = 0
            ],
            {k, Length[auxPairs]}
        ];

        eqList = Flatten[{
            List @@ (EqGeneral /. Thread[
                Keys[Lookup[Eqs, "costpluscurrents", <||>]] -> 0
            ]),
            If[Head[EqEntryIn] === And, List @@ EqEntryIn, If[EqEntryIn === True, {}, {EqEntryIn}]],
            If[Head[EqBalanceSplittingFlows] === And,
                List @@ EqBalanceSplittingFlows,
                If[EqBalanceSplittingFlows === True, {}, {EqBalanceSplittingFlows}]
            ],
            If[Head[EqValueAuxiliaryEdges] === And,
                List @@ EqValueAuxiliaryEdges,
                If[EqValueAuxiliaryEdges === True, {}, {EqValueAuxiliaryEdges}]
            ],
            KeyValueMap[Equal, RuleExitValues],
            KeyValueMap[Equal, RuleEntryValues],
            Equal @@@ Normal[RuleEntryOut],
            Equal @@@ Normal[RuleExitFlowsIn],
            Equal @@@ Normal[RuleBalanceGatheringFlows],
            (Equal[#, 0] & /@ Keys[zeroRules]),
            Module[{byVertex = GroupBy[auxPairs, Last]},
                Flatten @ KeyValueMap[
                    Function[{v, pairs},
                        If[
                            Length[pairs] > 1,
                            (u @@ # == u @@ pairs[[1]]) & /@ Rest[pairs],
                            {}
                        ]
                    ],
                    byVertex
                ]
            ]
        }];
        eqList = Cases[DeleteCases[eqList, True], _Equal];

        flowVars = Join[js, jts];
        allVars = Join[us, js, jts];
        ineqList = (# >= 0) & /@ flowVars;
        <|
            "Equalities" -> eqList,
            "Inequalities" -> ineqList,
            "FlowVariables" -> flowVars,
            "AllVariables" -> allVars
        |>
    ];

BuildLinearSystemFromEqualities[equalities_List, vars_List] :=
    Module[{exprs, placeholders, subExprs, coeffData, bVec, aMat},
        If[equalities === {} || vars === {},
            Return[$Failed, Module]
        ];
        exprs = equalities /. Equal[l_, r_] :> (l - r);
        placeholders = Table[Unique["lv"], {Length[vars]}];
        subExprs = exprs /. Thread[vars -> placeholders];
        coeffData = Quiet @ Check[
            CoefficientArrays[subExprs, placeholders],
            $Failed
        ];
        If[coeffData === $Failed || Length[coeffData] < 2,
            Return[$Failed, Module]
        ];
        bVec = Quiet @ Check[
            Developer`ToPackedArray @ N @ (-Normal[coeffData[[1]]]),
            $Failed
        ];
        aMat = Quiet @ Check[
            SparseArray[coeffData[[2]]],
            $Failed
        ];
        If[
            bVec === $Failed || aMat === $Failed || !VectorQ[bVec, NumericQ],
            $Failed,
            <|"A" -> aMat, "b" -> bVec|>
        ]
    ];

LinearSolveCandidate[aMat_SparseArray, bVec_List] :=
    Module[{dims, ls, x, ata, atb},
        dims = Dimensions[aMat];
        If[Length[dims] =!= 2 || dims[[2]] == 0,
            Return[$Failed, Module]
        ];
        If[dims[[1]] === dims[[2]],
            ls = Quiet @ Check[LinearSolve[N @ aMat], $Failed];
            If[ls === $Failed,
                Return[$Failed, Module]
            ];
            x = Quiet @ Check[N @ ls[N @ bVec], $Failed];
            If[VectorQ[x, NumericQ], Developer`ToPackedArray @ x, $Failed]
            ,
            ata = Quiet @ Check[N @ (Transpose[aMat].aMat), $Failed];
            atb = Quiet @ Check[N @ (Transpose[aMat].bVec), $Failed];
            If[ata === $Failed || atb === $Failed,
                Return[$Failed, Module]
            ];
            ls = Quiet @ Check[LinearSolve[ata], $Failed];
            If[ls === $Failed,
                x = Quiet @ Check[N @ LeastSquares[N @ Normal[aMat], N @ bVec], $Failed];
                Return[If[VectorQ[x, NumericQ], Developer`ToPackedArray @ x, $Failed], Module]
            ];
            x = Quiet @ Check[N @ ls[atb], $Failed];
            If[!VectorQ[x, NumericQ],
                x = Quiet @ Check[N @ LeastSquares[N @ Normal[aMat], N @ bVec], $Failed]
            ];
            If[VectorQ[x, NumericQ], Developer`ToPackedArray @ x, $Failed]
        ]
    ];

CriticalJFirstFailure[reason_String, details_Association : <||>] :=
    Failure["CriticalJFirstBackend", Join[<|"Reason" -> reason|>, details]];

CriticalJFirstFailureReason[result_] :=
    Which[
        MatchQ[result, Failure[_, _Association]],
            Lookup[result[[2]], "Reason", "JFirstFailed"],
        True,
            "JFirstFailed"
    ];

BuildCriticalJFirstObjectiveExpression[Eqs_Association, massVars_List] :=
    Module[{model, objective, objectiveVars},
        model = Quiet @ Check[BuildCriticalQuadraticObjective[Eqs], $Failed];
        If[!AssociationQ[model],
            Return[CriticalJFirstFailure["ObjectiveUnavailable"], Module]
        ];
        objective = Lookup[model, "Objective", Missing["NotAvailable"]];
        If[objective === Missing["NotAvailable"],
            Return[CriticalJFirstFailure["ObjectiveUnavailable"], Module]
        ];
        objectiveVars = Select[Variables[objective], MatchQ[#, _j] &];
        If[!SubsetQ[massVars, objectiveVars],
            Return[CriticalJFirstFailure["ObjectiveScopeMismatch"], Module]
        ];
        objective
    ];

CriticalConjunctionTerms[expr_] :=
    Which[
        TrueQ[expr], {},
        Head[expr] === And, List @@ expr,
        True, {expr}
    ];

RecoverCriticalFlowAssociation[Eqs_Association, flowBase_Association, tol_?NumericQ] :=
    Module[{entryInRules, flowRules, replacementRules, allFlowVars, resolve, resolved},
        entryInRules = Association @ Flatten[ToRules /@ Lookup[Eqs, "EqEntryIn", {}]];
        flowRules = Join[
            entryInRules,
            Lookup[Eqs, "RuleEntryOut", <||>],
            Lookup[Eqs, "RuleExitFlowsIn", <||>],
            Lookup[Eqs, "RuleBalanceGatheringFlows", <||>]
        ];
        replacementRules = Join[Normal[flowBase], Normal[flowRules]];
        resolve[expr_] := FixedPoint[
            Function[val, Quiet @ Check[val /. replacementRules, val]],
            expr,
            10
        ];
        allFlowVars = Join[Lookup[Eqs, "js", {}], Lookup[Eqs, "jts", {}]];
        resolved = Association @ Cases[
            Map[
                Function[{var},
                    Module[{val},
                        val = Quiet @ Check[N @ resolve[var], Missing["NotAvailable"]];
                        If[NumericQ[val],
                            var -> If[Abs[val] <= tol, 0., N[val]],
                            Nothing
                        ]
                    ]
                ],
                allFlowVars
            ],
            _Rule
        ];
        If[
            !AssociationQ[resolved] ||
            !And @@ (KeyExistsQ[resolved, #] & /@ allFlowVars),
            CriticalJFirstFailure["FlowReconstructionFailed"],
            resolved
        ]
    ];

SolveCriticalJFirstUtilities[
    Eqs_Association,
    flowAssoc_Association,
    uVars_List,
    tol_?NumericQ
] :=
    Module[{costpluscurrents, eqGeneralCritical, eqValueAuxiliaryEdges, ruleEntryValues,
      ruleExitValues, ineqSwitchingByVertex, switchingEqualities, uEqualities,
      uSystem, uCandidate, residual, utilityAssoc},
        If[uVars === {},
            Return[CriticalJFirstFailure["MissingUVars"], Module]
        ];
        costpluscurrents = Lookup[Eqs, "costpluscurrents", <||>];
        eqGeneralCritical =
            Lookup[Eqs, "EqGeneral", True] /.
                If[AssociationQ[costpluscurrents],
                    Expand /@ (costpluscurrents /. flowAssoc),
                    {}
                ];
        eqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", True];
        ruleEntryValues = Lookup[Eqs, "RuleEntryValues", <||>];
        ruleExitValues = Lookup[Eqs, "RuleExitValues", <||>];
        ineqSwitchingByVertex = Lookup[Eqs, "IneqSwitchingByVertex", True];
        switchingEqualities = DeleteDuplicatesBy[
            Cases[
                CriticalConjunctionTerms[ineqSwitchingByVertex],
                LessEqual[lhs_, rhs_] :> (lhs == rhs)
            ],
            Sort @* (List @@ # &)
        ];
        uEqualities = DeleteDuplicates @ Select[
            Cases[
                DeleteCases[
                    Flatten[{
                        CriticalConjunctionTerms[eqGeneralCritical],
                        CriticalConjunctionTerms[eqValueAuxiliaryEdges],
                        switchingEqualities,
                        KeyValueMap[Equal, ruleEntryValues],
                        KeyValueMap[Equal, ruleExitValues]
                    }] /. flowAssoc,
                    True
                ],
                _Equal
            ],
            !FreeQ[#, Alternatives @@ uVars] &
        ];
        If[uEqualities === {},
            Return[CriticalJFirstFailure["URecoveryFailed"], Module]
        ];
        uSystem = BuildLinearSystemFromEqualities[uEqualities, uVars];
        If[!AssociationQ[uSystem],
            Return[CriticalJFirstFailure["URecoveryFailed"], Module]
        ];
        uCandidate = LinearSolveCandidate[uSystem["A"], uSystem["b"]];
        If[!VectorQ[uCandidate, NumericQ] || Length[uCandidate] =!= Length[uVars],
            Return[CriticalJFirstFailure["URecoveryFailed"], Module]
        ];
        residual = Quiet @ Check[N @ Norm[uSystem["A"] . uCandidate - uSystem["b"], Infinity], Infinity];
        If[!NumericQ[residual] || residual > tol,
            Return[CriticalJFirstFailure["URecoveryFailed"], Module]
        ];
        utilityAssoc = AssociationThread[uVars, uCandidate];
        Association @ KeyValueMap[
            #1 -> If[Abs[#2] <= tol, 0., N[#2]] &,
            utilityAssoc
        ]
    ];

SolveCriticalJFirstBackend[Eqs_Association, numericState_Association] :=
    Module[{decoupling, topo, massConstraints, massVars, aEq, bEqRaw, bEq, tol = $CriticalSolverTolerance,
      nVars, massAssoc = <||>, xCandidate, candidateResidual, flowMin, rank,
      objectiveExpr, qpVars, qpRules, objectiveQP, constraintsQP,
      optimizerOrder, optimizer, rawSolution = $Failed, massVec,
      flowAssoc, uVars, utilityAssoc, solvedAssoc, expectedVars},
        decoupling = Lookup[numericState, "CriticalDecoupling", Missing["NotAvailable"]];
        If[!AssociationQ[decoupling],
            Return[CriticalJFirstFailure["MissingDecoupling"], Module]
        ];
        topo = Lookup[decoupling, "TopologicalOrder", <||>];
        If[!TrueQ[Lookup[topo, "IsDAG", False]],
            Return[CriticalJFirstFailure["NonDAG"], Module]
        ];
        massConstraints = Lookup[decoupling, "MassConstraints", <||>];
        massVars = Lookup[massConstraints, "Variables", {}];
        aEq = Lookup[massConstraints, "Aeq", Missing["NotAvailable"]];
        bEqRaw = Lookup[massConstraints, "beq", Missing["NotAvailable"]];
        If[Head[aEq] =!= SparseArray,
            aEq = SparseArray[aEq]
        ];
        bEq = Developer`ToPackedArray @ N @ Which[
            Head[bEqRaw] === SparseArray, Normal[bEqRaw],
            VectorQ[bEqRaw, NumericQ], bEqRaw,
            True, {}
        ];
        nVars = Length[massVars];
        If[
            nVars == 0 || Head[aEq] =!= SparseArray || !VectorQ[bEq, NumericQ] ||
            Length[bEq] =!= Dimensions[aEq][[1]] || Dimensions[aEq][[2]] =!= nVars,
            Return[CriticalJFirstFailure["MissingMassConstraints"], Module]
        ];

        rank = Quiet @ Check[MatrixRank[Normal[aEq]], -1];
        If[IntegerQ[rank] && rank === nVars,
            xCandidate = LinearSolveCandidate[aEq, bEq];
            If[VectorQ[xCandidate, NumericQ] && Length[xCandidate] === nVars,
                candidateResidual = Quiet @ Check[N @ Norm[aEq . xCandidate - bEq, Infinity], Infinity];
                flowMin = Min[xCandidate];
                If[
                    NumericQ[candidateResidual] && candidateResidual <= tol &&
                    NumericQ[flowMin] && flowMin >= -tol,
                    massAssoc = AssociationThread[massVars, xCandidate]
                ]
            ]
        ];

        If[!AssociationQ[massAssoc] || massAssoc === <||>,
            objectiveExpr = BuildCriticalJFirstObjectiveExpression[Eqs, massVars];
            qpVars = Table[Unique["qv"], {nVars}];
            qpRules = Thread[massVars -> qpVars];
            constraintsQP = Join[
                Thread[(aEq . qpVars) == bEq],
                Thread[qpVars >= 0]
            ];
            If[!MatchQ[objectiveExpr, Failure[_, _Association]],
                objectiveQP = objectiveExpr /. qpRules;
                optimizerOrder = DeleteDuplicates @ {
                    If[
                        NameQ["System`QuadraticOptimization"] &&
                        TrueQ @ Quiet @ Check[UseQuadraticCriticalBackendQ[Eqs], False],
                        QuadraticOptimization,
                        ConvexOptimization
                    ],
                    ConvexOptimization
                };
                Do[
                    rawSolution = Quiet @ Check[
                        optimizer[
                            objectiveQP,
                            constraintsQP,
                            qpVars
                        ],
                        $Failed
                    ];
                    If[MatchQ[rawSolution, {_Rule ..}],
                        Break[]
                    ],
                    {optimizer, optimizerOrder}
                ];
            ];
            If[!MatchQ[rawSolution, {_Rule ..}],
                rawSolution = Quiet @ Check[
                    LinearOptimization[
                        Total[qpVars],
                        constraintsQP,
                        qpVars
                    ],
                    $Failed
                ];
            ];
            If[!MatchQ[rawSolution, {_Rule ..}],
                Return[CriticalJFirstFailure["FlowInfeasible"], Module]
            ];
            massAssoc = AssociationThread[massVars, N @ (qpVars /. rawSolution)];
        ];

        massVec = Lookup[massAssoc, massVars, Missing["NotAvailable"]];
        If[!VectorQ[massVec, NumericQ] || Length[massVec] =!= nVars,
            Return[CriticalJFirstFailure["FlowInfeasible"], Module]
        ];
        candidateResidual = Quiet @ Check[N @ Norm[aEq . massVec - bEq, Infinity], Infinity];
        flowMin = Min[massVec];
        If[
            !NumericQ[candidateResidual] || candidateResidual > tol ||
            !NumericQ[flowMin] || flowMin < -tol,
            Return[CriticalJFirstFailure["FlowInfeasible"], Module]
        ];

        flowAssoc = RecoverCriticalFlowAssociation[Eqs, massAssoc, tol];
        If[MatchQ[flowAssoc, Failure[_, _Association]],
            Return[flowAssoc, Module]
        ];
        uVars = Lookup[decoupling, "UVars", Lookup[numericState, "UVars", {}]];
        utilityAssoc = SolveCriticalJFirstUtilities[Eqs, flowAssoc, uVars, tol];
        If[MatchQ[utilityAssoc, Failure[_, _Association]],
            Return[utilityAssoc, Module]
        ];
        solvedAssoc = Join[utilityAssoc, flowAssoc];
        expectedVars = Join[
            Lookup[Eqs, "us", {}],
            Lookup[Eqs, "js", {}],
            Lookup[Eqs, "jts", {}]
        ];
        If[!And @@ (KeyExistsQ[solvedAssoc, #] & /@ expectedVars),
            Return[CriticalJFirstFailure["IncompleteCandidate"], Module]
        ];
        Join[
            KeyTake[solvedAssoc, Lookup[Eqs, "us", {}]],
            KeyTake[solvedAssoc, Lookup[Eqs, "js", {}]],
            KeyTake[solvedAssoc, Lookup[Eqs, "jts", {}]]
        ]
    ];

SolveCriticalNumericBackendWithTelemetry[Eqs_Association] :=
    Module[{numericState, jFirstCandidate, jFirstReason, coupledCandidate,
      jFirstStatus, jFirstValidQ, forcedRejectQ},
        numericState = Lookup[Eqs, "NumericState", <||>];
        forcedRejectQ = TrueQ[Lookup[Eqs, "ForceNumericBackendValidationFailure", False]];
        jFirstCandidate = SolveCriticalJFirstBackend[Eqs, numericState];
        If[AssociationQ[jFirstCandidate],
            jFirstStatus = CheckFlowFeasibility[jFirstCandidate];
            jFirstValidQ =
                TrueQ[
                    !forcedRejectQ &&
                    jFirstStatus === "Feasible" &&
                    IsCriticalSolution[
                        Join[Eqs, <|"AssoCritical" -> jFirstCandidate, "Solution" -> jFirstCandidate|>],
                        "Tolerance" -> $CriticalSolverTolerance,
                        "Verbose" -> False
                    ]
                ];
            If[jFirstValidQ,
                Return[
                    <|
                        "Candidate" -> jFirstCandidate,
                        "Strategy" -> "JFirst",
                        "JFirstBackendUsed" -> True,
                        "JFirstBackendFallbackReason" -> None
                    |>,
                    Module
                ]
            ];
            jFirstReason = Which[
                forcedRejectQ, "ForcedValidationFailure",
                jFirstStatus =!= "Feasible", "FlowFeasibilityFailed",
                True, "CriticalValidationFailed"
            ],
            jFirstReason = CriticalJFirstFailureReason[jFirstCandidate]
        ];
        coupledCandidate = SolveCriticalNumericBackend[Eqs];
        <|
            "Candidate" -> coupledCandidate,
            "Strategy" -> "Coupled",
            "JFirstBackendUsed" -> False,
            "JFirstBackendFallbackReason" -> jFirstReason
        |>
    ];

SolveCriticalNumericBackend[Eqs_Association] :=
    Module[{constraints, equalities, flowVars, allVars,
        linearSystem, aMat, bVec, xCandidate, varIndex, flowIndices,
        linResidual, flowMin, tol = $CriticalSolverTolerance, assoc, lpResult, lpVals, solvedAssoc,
        nVars, cVec, aIneq, bIneq, aEq, bEq},
        constraints = BuildCriticalLinearConstraints[Eqs];
        If[constraints === $Failed || !AssociationQ[constraints],
            Return[$Failed, Module]
        ];
        equalities = Lookup[constraints, "Equalities", {}];
        flowVars = Lookup[constraints, "FlowVariables", {}];
        allVars = Lookup[constraints, "AllVariables", {}];
        If[equalities === {} || allVars === {},
            Return[$Failed, Module]
        ];
        varIndex = AssociationThread[allVars, Range[Length[allVars]]];
        flowIndices = Lookup[varIndex, flowVars, Missing["NotAvailable"]];
        If[MemberQ[flowIndices, _Missing] || flowIndices === {},
            Return[$Failed, Module]
        ];

        linearSystem = BuildLinearSystemFromEqualities[equalities, allVars];
        If[AssociationQ[linearSystem],
            aMat = linearSystem["A"];
            bVec = linearSystem["b"];
            xCandidate = LinearSolveCandidate[aMat, bVec];
            If[VectorQ[xCandidate, NumericQ] && Length[xCandidate] === Length[allVars],
                linResidual = Quiet @ Check[N @ Norm[aMat . xCandidate - bVec, Infinity], Infinity];
                flowMin = Min[xCandidate[[flowIndices]]];
                If[NumericQ[linResidual] && linResidual <= tol && NumericQ[flowMin] && flowMin >= -tol,
                    assoc = AssociationThread[allVars, xCandidate];
                    solvedAssoc = Association @ KeyValueMap[
                        #1 -> If[Abs[#2] <= tol, 0., N[#2]] &,
                        assoc
                    ];
                    Return[
                        Join[
                            KeyTake[solvedAssoc, Lookup[Eqs, "us", {}]],
                            KeyTake[solvedAssoc, Lookup[Eqs, "js", {}]],
                            KeyTake[solvedAssoc, Lookup[Eqs, "jts", {}]]
                        ],
                        Module
                    ]
                ]
            ]
            ,
            Return[$Failed, Module]
        ];

        nVars = Length[allVars];
        cVec = Developer`ToPackedArray @ ConstantArray[0., nVars];
        aEq = If[Head[aMat] === SparseArray, aMat, SparseArray[aMat]];
        bEq = Developer`ToPackedArray @ N @ (-bVec);
        aIneq = SparseArray[
            Thread[
                Transpose[{Range[Length[flowIndices]], flowIndices}] -> 1.
            ],
            {Length[flowIndices], nVars}
        ];
        bIneq = Developer`ToPackedArray @ ConstantArray[0., Length[flowIndices]];
        lpResult = Quiet @ Check[
            LinearOptimization[
                cVec,
                {aIneq, bIneq},
                {aEq, bEq},
                Tolerance -> tol
            ],
            $Failed
        ];
        lpVals = If[
            VectorQ[lpResult, NumericQ] && Length[lpResult] === nVars,
            Developer`ToPackedArray @ N @ lpResult,
            $Failed
        ];
        If[lpVals === $Failed,
            lpResult = Quiet @ Check[
                LinearOptimization[
                    cVec,
                    {aIneq, bIneq},
                    {aEq, bEq},
                    Method -> "Simplex",
                    Tolerance -> tol
                ],
                $Failed
            ];
            lpVals = If[
                VectorQ[lpResult, NumericQ] && Length[lpResult] === nVars,
                Developer`ToPackedArray @ N @ lpResult,
                $Failed
            ];
        ];
        If[
            lpVals === $Failed || !VectorQ[lpVals, NumericQ] || Length[lpVals] =!= nVars,
            Return[$Failed, Module]
        ];
        assoc = AssociationThread[allVars, Developer`ToPackedArray @ N @ lpVals];
        solvedAssoc = Association @ KeyValueMap[
            #1 -> If[Abs[#2] <= tol, 0., N[#2]] &,
            assoc
        ];
        Join[
            KeyTake[solvedAssoc, Lookup[Eqs, "us", {}]],
            KeyTake[solvedAssoc, Lookup[Eqs, "js", {}]],
            KeyTake[solvedAssoc, Lookup[Eqs, "jts", {}]]
        ]
    ];

(* --- CriticalCongestionSolver --- *)

Options[CriticalCongestionSolver] = {
    "ValidateSolution" -> True,
    "ValidationTolerance" -> $CriticalSolverTolerance,
    "ValidationVerbose" -> False
};

CriticalCongestionSolver[$Failed, ___] :=
    $Failed

CriticalCongestionSolver[Eqs_, OptionsPattern[]] :=
    Module[{PreEqs, js, AssoCritical, time, temp, status,
         resultKind, message, solution, comparisonData,
         validateQ, validationTolerance, validationVerboseQ, validationFailedQ,
         numericBackendRequestedQ, numericBackendEligibleQ,
         numericBackendUsedQ, numericBackendFallbackReason, numericBackendSolveTime,
         numericBackendStrategy, jFirstBackendUsedQ, jFirstBackendFallbackReason,
         numericCandidate, numericOutcome, forcedRejectQ, numericBackendTimeLimit},
        ClearSolveCache[];
        validateQ = TrueQ[OptionValue["ValidateSolution"]];
        validationTolerance = OptionValue["ValidationTolerance"];
        validationVerboseQ = TrueQ[OptionValue["ValidationVerbose"]];
        validationFailedQ = False;
        numericBackendUsedQ = False;
        numericBackendSolveTime = Missing["NotAvailable"];
        numericBackendStrategy = "Symbolic";
        jFirstBackendUsedQ = False;
        jFirstBackendFallbackReason = None;
        numericBackendRequestedQ = CriticalNumericBackendRequestedQ[Eqs];
        numericBackendEligibleQ = CriticalNumericBackendEligibleQ[Eqs];
        numericBackendFallbackReason = Which[
            !numericBackendRequestedQ, "DisabledByGuard",
            !numericBackendEligibleQ, "IneligibleInput",
            True, "NumericSolveFailed"
        ];
        forcedRejectQ = TrueQ[Lookup[Eqs, "ForceNumericBackendValidationFailure", False]];
        numericBackendTimeLimit = Lookup[Eqs, "CriticalNumericBackendTimeLimit", 30.];
        numericBackendTimeLimit = Which[
            numericBackendTimeLimit === Infinity || numericBackendTimeLimit === DirectedInfinity[1],
                Infinity,
            NumericQ[numericBackendTimeLimit] && numericBackendTimeLimit > 0,
                N[numericBackendTimeLimit],
            True,
                30.
        ];

        (* Try the numeric backend first when explicitly forced or globally enabled. *)
        If[numericBackendRequestedQ && numericBackendEligibleQ,
            MFGPrint["Attempting critical numeric backend..."];
            {numericBackendSolveTime, numericOutcome} =
                AbsoluteTiming[
                    TimeConstrained[
                        SolveCriticalNumericBackendWithTelemetry[Eqs],
                        numericBackendTimeLimit,
                        $TimedOut
                    ]
                ];
            If[numericOutcome === $TimedOut,
                numericBackendFallbackReason = "NumericBackendTimeout";
                numericBackendStrategy = "Symbolic";
                jFirstBackendUsedQ = False;
                jFirstBackendFallbackReason = "NumericBackendTimeout";
                MFGPrint[
                    "Numeric backend timed out after ",
                    numericBackendTimeLimit,
                    " seconds; falling back to symbolic path."
                ];
                numericCandidate = $Failed;
                ,
                numericCandidate = Lookup[numericOutcome, "Candidate", $Failed];
                numericBackendStrategy = Lookup[numericOutcome, "Strategy", "Coupled"];
                jFirstBackendUsedQ = TrueQ[Lookup[numericOutcome, "JFirstBackendUsed", False]];
                jFirstBackendFallbackReason = Lookup[
                    numericOutcome,
                    "JFirstBackendFallbackReason",
                    "NotAttempted"
                ];
            ];
            If[AssociationQ[numericCandidate] && numericCandidate =!= $Failed,
                status = CheckFlowFeasibility[numericCandidate];
                If[status === "Feasible",
                    numericBackendUsedQ = True;
                    If[KeyExistsQ[Eqs, "InitRules"],
                        PreEqs = Eqs,
                        {time, PreEqs} = AbsoluteTiming @ MFGPreprocessing[Eqs];
                        MFGPrint["Preprocessing (for downstream solvers) took ", time, " seconds."]
                    ];
                    If[validateQ,
                        If[forcedRejectQ || !TrueQ @ IsCriticalSolution[
                            Join[PreEqs, <|"AssoCritical" -> numericCandidate, "Solution" -> numericCandidate|>],
                            "Tolerance" -> validationTolerance,
                            "Verbose" -> validationVerboseQ
                        ],
                            numericBackendUsedQ = False;
                            numericBackendStrategy = "Symbolic";
                            If[
                                TrueQ[jFirstBackendUsedQ] &&
                                (jFirstBackendFallbackReason === None || jFirstBackendFallbackReason === "NotAttempted"),
                                jFirstBackendFallbackReason =
                                    If[forcedRejectQ, "ForcedValidationFailure", "CriticalValidationFailed"]
                            ];
                            jFirstBackendUsedQ = False;
                            numericBackendFallbackReason =
                                If[forcedRejectQ, "ForcedValidationFailure", "CriticalValidationFailed"]
                            ,
                            solution = numericCandidate;
                            comparisonData = BuildSolverComparisonData[PreEqs, solution];
                            Return[
                                Join[
                                    PreEqs,
                                    MakeSolverResult[
                                        "CriticalCongestion",
                                        "Success",
                                        status,
                                        None,
                                        solution,
                                        Join[
                                            comparisonData,
                                            <|
                                                "AssoCritical" -> numericCandidate,
                                                "Status" -> status,
                                                "NumericBackendUsed" -> True,
                                                "NumericBackendFallbackReason" -> None,
                                                "NumericBackendSolveTime" -> numericBackendSolveTime,
                                                "NumericBackendStrategy" -> numericBackendStrategy,
                                                "JFirstBackendUsed" -> jFirstBackendUsedQ,
                                                "JFirstBackendFallbackReason" -> jFirstBackendFallbackReason
                                            |>
                                        ]
                                    ]
                                ],
                                Module
                            ]
                        ],
                        solution = numericCandidate;
                        comparisonData = BuildSolverComparisonData[PreEqs, solution];
                        Return[
                            Join[
                                PreEqs,
                                MakeSolverResult[
                                    "CriticalCongestion",
                                    "Success",
                                    status,
                                    None,
                                    solution,
                                    Join[
                                        comparisonData,
                                        <|
                                            "AssoCritical" -> numericCandidate,
                                            "Status" -> status,
                                            "NumericBackendUsed" -> True,
                                            "NumericBackendFallbackReason" -> None,
                                            "NumericBackendSolveTime" -> numericBackendSolveTime,
                                            "NumericBackendStrategy" -> numericBackendStrategy,
                                            "JFirstBackendUsed" -> jFirstBackendUsedQ,
                                            "JFirstBackendFallbackReason" -> jFirstBackendFallbackReason
                                        |>
                                    ]
                                ]
                            ],
                            Module
                        ]
                    ],
                    numericBackendFallbackReason = "FlowFeasibilityFailed"
                ],
                If[numericBackendFallbackReason =!= "NumericBackendTimeout",
                    numericBackendFallbackReason = "NumericSolveFailed"
                ]
            ]
        ];

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
                    If[validateQ,
                        If[!TrueQ @ IsCriticalSolution[
                            Join[PreEqs, <|"AssoCritical" -> AssoCritical, "Solution" -> AssoCritical|>],
                            "Tolerance" -> validationTolerance,
                            "Verbose" -> validationVerboseQ
                        ],
                            MFGPrint["Direct solver validation failed, falling back to symbolic solver."]
                            ,
                            solution = AssoCritical;
                            comparisonData = BuildSolverComparisonData[PreEqs, solution];
                            Return[
                                Join[PreEqs, MakeSolverResult["CriticalCongestion", "Success",
                                    status, None, solution,
                                    Join[
                                        comparisonData,
                                        <|
                                            "AssoCritical" -> AssoCritical,
                                            "Status" -> status,
                                            "NumericBackendUsed" -> numericBackendUsedQ,
                                            "NumericBackendFallbackReason" -> numericBackendFallbackReason,
                                            "NumericBackendSolveTime" -> numericBackendSolveTime,
                                            "NumericBackendStrategy" -> numericBackendStrategy,
                                            "JFirstBackendUsed" -> jFirstBackendUsedQ,
                                            "JFirstBackendFallbackReason" -> jFirstBackendFallbackReason
                                        |>
                                    ]]],
                                Module
                            ]
                        ],
                        solution = AssoCritical;
                        comparisonData = BuildSolverComparisonData[PreEqs, solution];
                        Return[
                            Join[PreEqs, MakeSolverResult["CriticalCongestion", "Success",
                                status, None, solution,
                                Join[
                                    comparisonData,
                                    <|
                                        "AssoCritical" -> AssoCritical,
                                        "Status" -> status,
                                        "NumericBackendUsed" -> numericBackendUsedQ,
                                        "NumericBackendFallbackReason" -> numericBackendFallbackReason,
                                        "NumericBackendSolveTime" -> numericBackendSolveTime,
                                        "NumericBackendStrategy" -> numericBackendStrategy,
                                        "JFirstBackendUsed" -> jFirstBackendUsedQ,
                                        "JFirstBackendFallbackReason" -> jFirstBackendFallbackReason
                                    |>
                                ]]],
                            Module
                        ]
                    ],
                    MFGPrint["Direct solver produced infeasible result, falling back to symbolic solver."]
                ],
                MFGPrint["Direct solver failed, falling back to symbolic solver."]
            ]
        ];
        (* Standard symbolic solver path *)
        If[!(AssociationQ[PreEqs] && KeyExistsQ[PreEqs, "InitRules"]),
            If[KeyExistsQ[Eqs, "InitRules"],
                PreEqs = Eqs
                ,
                temp = MFGPrintTemporary["Preprocessing..."];
                {time, PreEqs} = AbsoluteTiming @ MFGPreprocessing[Eqs];
                NotebookDelete[temp];
                MFGPrint["Preprocessing took ", time, " seconds to terminate."
                    ];
            ]
        ];
        js = Lookup[PreEqs, "js", $Failed];
        AssoCritical = MFGSystemSolver[PreEqs][AssociationThread[js,
            0 js]];
        status = CheckFlowFeasibility[AssoCritical];
        If[validateQ && AssociationQ[AssoCritical] && status === "Feasible",
            If[!TrueQ @ IsCriticalSolution[
                Join[PreEqs, <|"AssoCritical" -> AssoCritical, "Solution" -> AssoCritical|>],
                "Tolerance" -> validationTolerance,
                "Verbose" -> validationVerboseQ
            ],
                validationFailedQ = True;
                AssoCritical = Null;
                status = "Infeasible";
            ]
        ];
        resultKind =
            If[AssoCritical === Null,
                "Failure"
                ,
                "Success"
            ];
        message =
            If[AssoCritical === Null,
                If[validationFailedQ, "InvalidCriticalSolution", "NoSolution"]
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
             status, message, solution,
             Join[
                 comparisonData,
                 <|
                    "AssoCritical" -> AssoCritical,
                    "Status" -> status,
                    "NumericBackendUsed" -> numericBackendUsedQ,
                    "NumericBackendFallbackReason" -> numericBackendFallbackReason,
                    "NumericBackendSolveTime" -> numericBackendSolveTime,
                    "NumericBackendStrategy" -> numericBackendStrategy,
                    "JFirstBackendUsed" -> jFirstBackendUsedQ,
                    "JFirstBackendFallbackReason" -> jFirstBackendFallbackReason
                 |>
             ]]]
    ];

Options[IsCriticalSolution] = {
    "Tolerance" -> $CriticalSolverTolerance,
    "Verbose" -> False,
    "ReturnReport" -> False
};

NumericRelationSatisfiedQ[expr_, tol_?NumericQ] :=
    Module[{delta, pair, held},
        held = HoldComplete[expr];
        Which[
            MatchQ[held, HoldComplete[Equal[_, _]] | HoldComplete[Inactive[Equal][_, _]]],
                pair = Replace[held, {
                    HoldComplete[Equal[l_, r_]] :> {l, r},
                    HoldComplete[Inactive[Equal][l_, r_]] :> {l, r}
                }];
                delta = Quiet @ Check[N[pair[[1]] - pair[[2]]], $Failed];
                NumericQ[delta] && Abs[delta] <= tol
            ,
            MatchQ[held, HoldComplete[LessEqual[_, _]] | HoldComplete[Inactive[LessEqual][_, _]]],
                pair = Replace[held, {
                    HoldComplete[LessEqual[l_, r_]] :> {l, r},
                    HoldComplete[Inactive[LessEqual][l_, r_]] :> {l, r}
                }];
                delta = Quiet @ Check[N[pair[[1]] - pair[[2]]], $Failed];
                NumericQ[delta] && delta <= tol
            ,
            MatchQ[held, HoldComplete[GreaterEqual[_, _]] | HoldComplete[Inactive[GreaterEqual][_, _]]],
                pair = Replace[held, {
                    HoldComplete[GreaterEqual[l_, r_]] :> {l, r},
                    HoldComplete[Inactive[GreaterEqual][l_, r_]] :> {l, r}
                }];
                delta = Quiet @ Check[N[pair[[2]] - pair[[1]]], $Failed];
                NumericQ[delta] && delta <= tol
            ,
            MatchQ[held, HoldComplete[Less[_, _]] | HoldComplete[Inactive[Less][_, _]]],
                pair = Replace[held, {
                    HoldComplete[Less[l_, r_]] :> {l, r},
                    HoldComplete[Inactive[Less][l_, r_]] :> {l, r}
                }];
                delta = Quiet @ Check[N[pair[[1]] - pair[[2]]], $Failed];
                NumericQ[delta] && delta < -tol
            ,
            MatchQ[held, HoldComplete[Greater[_, _]] | HoldComplete[Inactive[Greater][_, _]]],
                pair = Replace[held, {
                    HoldComplete[Greater[l_, r_]] :> {l, r},
                    HoldComplete[Inactive[Greater][l_, r_]] :> {l, r}
                }];
                delta = Quiet @ Check[N[pair[[2]] - pair[[1]]], $Failed];
                NumericQ[delta] && delta < -tol
            ,
            MatchQ[held, HoldComplete[Unequal[_, _]] | HoldComplete[Inactive[Unequal][_, _]]],
                pair = Replace[held, {
                    HoldComplete[Unequal[l_, r_]] :> {l, r},
                    HoldComplete[Inactive[Unequal][l_, r_]] :> {l, r}
                }];
                delta = Quiet @ Check[N[pair[[1]] - pair[[2]]], $Failed];
                NumericQ[delta] && Abs[delta] > tol
            ,
            True,
                False
        ]
    ];

LogicalSatisfiedQ[expr_, rules_Association, tol_?NumericQ] :=
    Module[{heldExpr, evaluated},
        Which[
            expr === True, True,
            expr === False, False,
            True,
                heldExpr =
                    expr /. HoldPattern[(h:Equal | LessEqual | GreaterEqual | Less | Greater | Unequal)[lhs_, rhs_]] :>
                        Inactive[h][lhs, rhs];
                evaluated = Quiet @ Check[
                    heldExpr /. rules /. {
                        HoldPattern[Inactive[Equal][lhs_, rhs_]] :>
                            NumericRelationSatisfiedQ[Inactive[Equal][lhs, rhs], tol],
                        HoldPattern[Inactive[LessEqual][lhs_, rhs_]] :>
                            NumericRelationSatisfiedQ[Inactive[LessEqual][lhs, rhs], tol],
                        HoldPattern[Inactive[GreaterEqual][lhs_, rhs_]] :>
                            NumericRelationSatisfiedQ[Inactive[GreaterEqual][lhs, rhs], tol],
                        HoldPattern[Inactive[Less][lhs_, rhs_]] :>
                            NumericRelationSatisfiedQ[Inactive[Less][lhs, rhs], tol],
                        HoldPattern[Inactive[Greater][lhs_, rhs_]] :>
                            NumericRelationSatisfiedQ[Inactive[Greater][lhs, rhs], tol],
                        HoldPattern[Inactive[Unequal][lhs_, rhs_]] :>
                            NumericRelationSatisfiedQ[Inactive[Unequal][lhs, rhs], tol]
                    },
                    $Failed
                ];
                If[evaluated === $Failed,
                    False,
                    TrueQ[Simplify[evaluated]]
                ]
        ]
    ];

RulesConjunction[rules_] :=
    Module[{pairs},
        pairs = Which[
            AssociationQ[rules], Normal[rules],
            ListQ[rules] && VectorQ[rules, MatchQ[#, _Rule]&], rules,
            True, {}
        ];
        If[pairs === {},
            True,
            And @@ (Equal @@@ (List @@@ pairs))
        ]
    ];

IsCriticalSolution[Eqs_Association, OptionsPattern[]] :=
    Module[{tol, verboseQ, returnReportQ, solution, blocks, blockResults,
         residual, residualPass, allBlocksPass, overall, report, flowApprox,
         costpluscurrents, eqGeneralCritical, eqResidualPairs, eqResidualDeltas,
         requiredKeys, missingKeys},
        tol = OptionValue["Tolerance"];
        verboseQ = TrueQ[OptionValue["Verbose"]];
        returnReportQ = TrueQ[OptionValue["ReturnReport"]];

        requiredKeys = {
            "EqEntryIn",
            "EqValueAuxiliaryEdges",
            "AltOptCond",
            "EqBalanceSplittingFlows",
            "AltFlows",
            "AltTransitionFlows",
            "IneqJs",
            "IneqJts",
            "IneqSwitchingByVertex",
            "RuleEntryOut",
            "RuleExitFlowsIn",
            "RuleEntryValues",
            "RuleExitValues",
            "RuleBalanceGatheringFlows",
            "EqGeneral",
            "js",
            "costpluscurrents"
        };
        missingKeys = Select[requiredKeys, !KeyExistsQ[Eqs, #]&];
        If[missingKeys =!= {},
            report = <|
                "Valid" -> False,
                "Reason" -> "MissingRequiredFields",
                "MissingKeys" -> missingKeys,
                "BlockChecks" -> <||>,
                "EquationResidual" -> Missing["NotAvailable"],
                "EquationResidualPass" -> False,
                "Tolerance" -> tol
            |>;
            If[verboseQ,
                Print["IsCriticalSolution: missing required fields: ", missingKeys]
            ];
            Return[If[returnReportQ, report, False], Module]
        ];

        solution = Which[
            AssociationQ[Lookup[Eqs, "AssoCritical", Missing["NotAvailable"]]],
                Lookup[Eqs, "AssoCritical"],
            AssociationQ[Lookup[Eqs, "Solution", Missing["NotAvailable"]]],
                Lookup[Eqs, "Solution"],
            True,
                Missing["NotAvailable"]
        ];

        If[!AssociationQ[solution],
            report = <|
                "Valid" -> False,
                "Reason" -> "MissingAssoCritical",
                "MissingKeys" -> {},
                "BlockChecks" -> <||>,
                "EquationResidual" -> Missing["NotAvailable"],
                "EquationResidualPass" -> False,
                "Tolerance" -> tol
            |>;
            If[verboseQ, Print["IsCriticalSolution: missing association in \"AssoCritical\" or \"Solution\"."]];
            Return[If[returnReportQ, report, False], Module]
        ];

        flowApprox = Lookup[Eqs, "js", {}];
        flowApprox = AssociationThread[flowApprox, 0 flowApprox];
        costpluscurrents = Lookup[Eqs, "costpluscurrents", Missing["NotAvailable"]];
        eqGeneralCritical =
            Lookup[Eqs, "EqGeneral", True] /.
                If[AssociationQ[costpluscurrents],
                    Expand /@ (costpluscurrents /. flowApprox),
                    {}
                ];

        blocks = <|
            "EqEntryIn" -> (And @@ Lookup[Eqs, "EqEntryIn", True]),
            "EqValueAuxiliaryEdges" -> Lookup[Eqs, "EqValueAuxiliaryEdges", True],
            "AltOptCond" -> Lookup[Eqs, "AltOptCond", True],
            "EqBalanceSplittingFlows" -> Lookup[Eqs, "EqBalanceSplittingFlows", True],
            "AltFlows" -> Lookup[Eqs, "AltFlows", True],
            "AltTransitionFlows" -> Lookup[Eqs, "AltTransitionFlows", True],
            "IneqJs" -> Lookup[Eqs, "IneqJs", True],
            "IneqJts" -> Lookup[Eqs, "IneqJts", True],
            "IneqSwitchingByVertex" -> (And @@ Lookup[Eqs, "IneqSwitchingByVertex", True]),
            "RuleEntryOut" -> RulesConjunction[Lookup[Eqs, "RuleEntryOut", <||>]],
            "RuleExitFlowsIn" -> RulesConjunction[Lookup[Eqs, "RuleExitFlowsIn", <||>]],
            "RuleEntryValues" -> RulesConjunction[Lookup[Eqs, "RuleEntryValues", <||>]],
            "RuleExitValues" -> RulesConjunction[Lookup[Eqs, "RuleExitValues", <||>]],
            "RuleBalanceGatheringFlows" -> RulesConjunction[Lookup[Eqs, "RuleBalanceGatheringFlows", <||>]],
            "EqGeneralCritical" -> eqGeneralCritical
        |>;

        blockResults = Association @ KeyValueMap[
            #1 -> LogicalSatisfiedQ[#2, solution, tol] &,
            blocks
        ];
        allBlocksPass = And @@ Values[blockResults];

        eqResidualPairs =
            Cases[
                (eqGeneralCritical /. HoldPattern[Equal[lhs_, rhs_]] :> Inactive[Equal][lhs, rhs]) /. solution,
                HoldPattern[Inactive[Equal][lhs_, rhs_]] :> {lhs, rhs},
                Infinity
            ];
        eqResidualDeltas =
            Quiet @ Check[
                N[(#[[1]] - #[[2]])& /@ eqResidualPairs],
                $Failed
            ];
        residual =
            If[eqResidualPairs === {},
                Missing["NotAvailable"],
                If[eqResidualDeltas === $Failed || !VectorQ[eqResidualDeltas, NumericQ],
                    Missing["ComputeError"],
                    Max[Abs[eqResidualDeltas]]
                ]
            ];
        residualPass =
            Which[
                NumericQ[residual], residual <= tol,
                residual === Missing["NotAvailable"], TrueQ[Lookup[blockResults, "EqGeneralCritical", False]],
                True, False
            ];

        overall = allBlocksPass && residualPass;

        If[verboseQ,
            Print["IsCriticalSolution: all blocks pass = ", allBlocksPass];
            If[!allBlocksPass,
                Print["  Failed blocks: ", Keys @ Select[blockResults, !TrueQ[#]&]]
            ];
            Print["IsCriticalSolution: equation residual = ", residual, " (tol=", tol, ")"];
            Print["IsCriticalSolution: equation residual pass = ", residualPass];
        ];

        report = <|
            "Valid" -> overall,
            "Reason" -> If[overall, None, "ConstraintViolation"],
            "MissingKeys" -> {},
            "BlockChecks" -> blockResults,
            "EquationResidual" -> residual,
            "EquationResidualPass" -> residualPass,
            "Tolerance" -> tol
        |>;

        If[returnReportQ, report, overall]
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
