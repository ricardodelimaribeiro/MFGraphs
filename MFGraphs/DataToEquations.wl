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



MFGSystemSolver::usage = "MFGSystemSolver[Eqs][edgeEquations] returns an Association \
<|\"Solution\" -> assoc_or_Null, \"UnresolvedConstraints\" -> expr_or_None|>. \
\"Solution\" is Null on failure; \"UnresolvedConstraints\" is non-None when flows were determined \
but u-variables remain symbolically underdetermined.";

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
         blockedIncomingPairs, activeAuxTriples,
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
                Message[DataToEquations::switchingcosts]
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
        (* Transition flow constraints: for each vertex pair {a,b}, at least one
           transition through {a,b} must have zero flow (complementarity).

           When switching costs are CONSISTENT (satisfy triangle inequality):
           Use original method: AltFlowOp on ordered vertex pairs only.

           When switching costs are INCONSISTENT (violate triangle inequality):
           Use disjunctive constraint strategy: For each target vertex v, for each
           pair of incoming transitions (t1, t2), require j(t1)==0 OR j(t2)==0.
           Groups auxTriples by intermediate vertex (midpoint), then for each target,
           constructs disjunctions among all paths reaching that target.
           This weaker constraint allows solver to find solutions when triangle
           inequality fails, while still enforcing some complementarity. *)
        AltTransitionFlows =
            If[consistentCosts === False,
                And @@ DeleteDuplicates[
                    Sort /@ Flatten @ KeyValueMap[
                        Function[{k, trips},
                            Module[{bySource, byTarget},
                                bySource = GroupBy[trips, First];
                                byTarget = GroupBy[trips, Last];
                                KeyValueMap[
                                    Function[{v, t1s},
                                        Table[
                                            j @@ t1 == 0 || j @@ t2 == 0,
                                            {t1, t1s},
                                            {t2, Lookup[bySource, v, {}]}
                                        ]
                                    ],
                                    byTarget
                                ]
                            ]
                        ],
                        GroupBy[auxTriples, #[[2]] &]
                    ]
                ],
                And @@ (AltFlowOp[j] /@
                    Select[auxTriples, OrderedQ[{First[#], Last[#]}]&])
            ];
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
        blockedIncomingPairs = DeleteDuplicates @ Join[outAuxEntryPairs, inAuxExitPairs];
        activeAuxTriples = Select[
            auxTriples,
            !MemberQ[blockedIncomingPairs, #[[{1, 2}]]] &
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
        IneqSwitchingByVertex =
            IneqSwitch[u, SwitchingCosts] @@@
                Select[
                    TransitionsAt[auxiliaryGraph, #],
                    !MemberQ[blockedIncomingPairs, #[[{1, 2}]]] &
                ] & /@ verticesList;
        IneqSwitchingByVertex = And @@@ IneqSwitchingByVertex;
        AltOptCond = And @@ AltSwitch[j, u, SwitchingCosts] @@@ activeAuxTriples
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
        RuleEntryValues, time},
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
                Module[{items = IneqSwitchingByVertex /. InitRules, ineqs, pairs,
                        reciprocalPairs, eqs, canonicalPair, reciprocalKeys,
                        remainingIneqs, combinedConditions},
                    ineqs = Flatten[If[Head[#] === And, List @@ #, {#}]& /@ items];
                    pairs = Cases[ineqs, LessEqual[a_, b_] :> {a, b}];
                    reciprocalPairs = Select[pairs, MemberQ[pairs, Reverse[#]] &];
                    canonicalPair[pair_] := SortBy[pair, ToString[InputForm[#]] &];
                    eqs = DeleteDuplicatesBy[Equal @@@ reciprocalPairs, canonicalPair @* (List @@ # &)];
                    reciprocalKeys = DeleteDuplicates[canonicalPair /@ reciprocalPairs];
                    remainingIneqs = Select[
                        ineqs,
                        Not @ MatchQ[
                            #,
                            LessEqual[a_, b_] /; MemberQ[reciprocalKeys, canonicalPair[{a, b}]]
                        ] &
                    ];
                    combinedConditions = Join[eqs, remainingIneqs];
                    If[combinedConditions === {}, True, And @@ combinedConditions]
                ]
            ];
            MFGPrint["Switching costs (zero): extracted equalities in ",
                time, " seconds."];
            ,
            With[{items = IneqSwitchingByVertex /. InitRules},
                {time, IneqSwitchingByVertex} = AbsoluteTiming[And @@ MFGParallelMap[
                    Simplify, items]]
            ];
            MFGPrint["Switching costs simplified in ", time, " seconds."];
        ];
        EqBalanceSplittingFlows = EqBalanceSplittingFlows /. InitRules
            ;
        EqValueAuxiliaryEdges = EqValueAuxiliaryEdges /. InitRules;
        {time, Rules} = AbsoluteTiming[First @ Solve[EqBalanceSplittingFlows
             && EqValueAuxiliaryEdges] // Quiet];
        InitRules = Join[InitRules /. Rules, Association @ Rules];
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
        (* Disabled: GraphDistance heuristic makes wrong zero-flow assumptions
           for heterogeneous-exit-cost cases (e.g. chain with two exits). *)
        (* {time, NewSystem} = AbsoluteTiming[
            ResolveOrByGraphDistance[Eqs, NewSystem]];
        MFGPrint["Graph-distance Or-resolution took ", time, " seconds."]; *)
        {time, {NewSystem, InitRules}} = AbsoluteTiming @ TripleClean[
            {NewSystem, InitRules}];
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

(* --- CriticalCongestionSolver --- *)

(* BuildPrunedSystem: given a numeric candidate solution and the symbolic system {EE, NN, OR},
   retain only the Or-branches satisfied by the candidate. Falls back to the full system
   if no branches are satisfied. *)
BuildPrunedSystem[system_, candidateSolution_Association, tol_: 10^-6] :=
    Module[{ee, nn, orBranches, prunedOr},
        {ee, nn, orBranches} = system;
        If[orBranches === True,
            prunedOr = True,
            prunedOr = If[Head[orBranches] === List, Or @@ orBranches, orBranches];
            prunedOr = prunedOr /. HoldPattern[Or[args__]] :>
                Module[{satisfied, lst = {args}},
                    satisfied = Select[lst, TrueQ[LogicalSatisfiedQ[#, candidateSolution, tol]] &];
                    If[satisfied === {}, Or[args], If[Length[satisfied] === 1, satisfied[[1]], Or @@ satisfied]]
                ];
            If[Head[orBranches] === List,
                prunedOr = If[Head[prunedOr] === Or, List @@ prunedOr, {prunedOr}]
            ]
        ];
        {ee, nn, prunedOr}
    ];

(* BuildOraclePrunedSystem: Translate OracleState payload from ClassifyAndCheckStability
   into a partially pruned {EE, NN, OR} triple. Uses two complementary strategies:
   1. LogicalSatisfiedQ-based OR pruning (existing BuildPrunedSystem machinery)
   2. Injection of explicit zero equalities into EE for confirmed-inactive variables

   PRECISION WARNING: The EE injection of j[e]==0 relies on floating-point thresholds
   (10^-5 for flow). If true equilibrium requires j[e] slightly below this threshold,
   over-constraining EE may cause symbolic solve to fail. See project memory
   "Phase 5/6 Technical Debt" for mitigation strategy (safety bands).

   Safe because PrunedZeroFlows contains only variables classified inactive by rigorous
   hysteresis thresholds. Ambiguous variables are untouched. *)
BuildOraclePrunedSystem[system_, flowAssoc_Association, oracleState_Association, tol_: 10^-6] :=
    Module[{prunedSystem, ee, nn, or, inactiveVars, newEE},
        (* Step 1: filter OR branches using existing LogicalSatisfiedQ machinery *)
        prunedSystem = BuildPrunedSystem[system, flowAssoc, tol];
        {ee, nn, or} = prunedSystem;

        (* Step 2: inject explicit zero equalities for confirmed-inactive variables into EE.
           This lets TripleClean propagate substitutions before DNFSolveStep runs,
           compounding the OR-branch reduction from step 1. *)
        inactiveVars = Lookup[oracleState, "PrunedZeroFlows", {}];
        newEE = If[ListQ[inactiveVars] && inactiveVars =!= {},
            ee && (And @@ (# == 0 & /@ inactiveVars)),
            ee
        ];
        {newEE, nn, or}
    ];

(* BuildCriticalResult: single point for assembling the standardized result envelope.
   Eliminates duplication across numeric, direct, and symbolic solver paths. *)
(* --- MFGSystemSolver --- *)

MFGSystemSolver[Eqs_][approxJs_] :=
    Module[{NewSystem, InitRules, pickOne, pickOneFlowOnly, vars, System, Ncpc, costpluscurrents,
         us, js, jts, jjtsR, time, ineqsByTransition, uResidual,
         ineqsWithoutTransition = {}},
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
            Return[<|"Solution" -> Null, "UnresolvedConstraints" -> None|>, Module]
        ];
(* Retrieve some equalities from the inequalities: group by transition flow 
    
    
    
    *)
        If[NewSystem[[2]] === True,
            ineqsByTransition = ConstantArray[True, Length[jts]];
            ineqsWithoutTransition = {};
            time = 0.
            ,
            {time, ineqsByTransition} =
                AbsoluteTiming[
                    Module[{ineqList, index, scopedIneqs, transitionPattern},
                        ineqList = If[Head[NewSystem[[2]]] === And,
                            List @@ NewSystem[[2]], {NewSystem[[2]]}];
                        transitionPattern = If[jts === {}, None, Alternatives @@ jts];
                        ineqsWithoutTransition =
                            If[transitionPattern === None,
                                ineqList,
                                Select[ineqList, FreeQ[#, transitionPattern] &]
                            ];
                        scopedIneqs =
                            If[ineqsWithoutTransition === {},
                                ineqList,
                                Complement[ineqList, ineqsWithoutTransition]
                            ];
                        index = Association[# -> {} & /@ jts];
                        Do[
                            Do[
                                If[!FreeQ[ineq, jt],
                                    index[jt] = Append[index[jt], ineq]],
                                {jt, jts}
                            ],
                            {ineq, scopedIneqs}
                        ];
                        (If[# === {}, True, And @@ #]&) /@ Values[index]
                    ]
                ]
        ];
        MFGPrint["MFGSS: Selecting inequalities by transition flow took ",
             time, " seconds. ", Length[ineqsByTransition]];
        (* Pre-substitute known rules: many inequalities become True immediately *)
        ineqsByTransition = ineqsByTransition /. InitRules;
        ineqsByTransition = Replace[ineqsByTransition,
            expr_ /; expr =!= True :> Expand[expr], {1}];
        {time, ineqsByTransition} = AbsoluteTiming[DeduplicateByComplexity[
            MFGParallelMap[Simplify,
                Select[ineqsByTransition, # =!= True &]]]];
        {time, ineqsWithoutTransition} =
            AbsoluteTiming[
                DeduplicateByComplexity[
                    MFGParallelMap[Simplify,
                        Select[ineqsWithoutTransition, # =!= True &]]
                ]
            ];
        (* Restore True entries for any that were already resolved *)
        MFGPrint["MFGSS: Simplifying inequalities by transition flow took ",
             time, " seconds. "];
        NewSystem[[2]] = (And @@ Join[ineqsByTransition, ineqsWithoutTransition]);
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
                Return[<|"Solution" -> Null, "UnresolvedConstraints" -> None|>, Module]
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
                        (* Both full-variable attempts failed. Try flows-only: WL existentially
                           quantifies over u vars, finding j values for which some u satisfies System. *)
                        pickOneFlowOnly = If[jjtsR =!= {},
                            FindInstance[System && And @@ ((# >= 0)& /@ jjtsR), jjtsR, Reals],
                            {}
                        ];
                        If[pickOneFlowOnly =!= {},
                            pickOneFlowOnly = Association @ First @ pickOneFlowOnly;
                            MFGPrint["MFGSS: Flow-only solution found (u underdetermined): ", pickOneFlowOnly];
                            InitRules = Expand /@ Join[InitRules /. pickOneFlowOnly, pickOneFlowOnly];
                            uResidual = Simplify[System /. pickOneFlowOnly];
                            (* Resolve transitive chains early so the result is usable downstream *)
                            InitRules = Expand /@ FixedPoint[Function[r, ReplaceAll[r] /@ r], InitRules, 10];
                            InitRules = Join[KeyTake[InitRules, us], KeyTake[InitRules, js], KeyTake[InitRules, jts]];
                            Return[<|"Solution" -> InitRules, "UnresolvedConstraints" -> uResidual|>, Module]
                            ,
                            MFGPrint["MFGSS: No feasible solution found"];
                            Return[<|"Solution" -> Null, "UnresolvedConstraints" -> None|>, Module]
                        ]
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
        <|"Solution" -> InitRules, "UnresolvedConstraints" -> None|>
    ];

(* --- DNFSolveStep --- *)

DNFSolveStep[{{EE_?BooleanQ, NN_?BooleanQ, OR_?BooleanQ}, rules_}] :=
    {{EE, NN, OR}, rules}

DNFSolveStep[{{EE_, NN_, OO_}, rules_}] :=
    Module[{NewSystem, newrules, sorted = True, time},
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
        {time, NewSystem} = AbsoluteTiming[DNFReduce[And @@ Most[NewSystem
            ], sorted]];
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
