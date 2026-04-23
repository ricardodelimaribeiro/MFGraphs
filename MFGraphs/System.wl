(* Wolfram Language package *)
(*
   mfgSystem: typed system kernel for MFGraphs.

   Encapsulates the structural equations, inequalities, and substitution rules
   that represent a stationary Mean Field Game on a network. This kernel
   separates symbolic formulation from the solver lifecycle.

   Lifecycle:
     makeSystem[s, unk]  →  builds system association + wraps in mfgSystem[...]
     makeSystem[s]       →  calls makeUnknowns[s] and delegates
     mfgSystemQ[x]       →  True iff x is a well-formed typed mfgSystem
     SystemData[sys]     →  returns the underlying equation association
*)

(* --- Public API declarations --- *)

mfgSystem::usage =
"mfgSystem[assoc] is the typed head for an MFG structural equation system. Use \
makeSystem to construct one; use SystemData to access keys.";

mfgSystem::switchingcosts = "Switching costs are inconsistent.";

mfgSystemQ::usage =
"mfgSystemQ[x] returns True if x is a typed mfgSystem[assoc_Association] object, \
False otherwise.";

makeSystem::usage =
"makeSystem[s_scenario, unk_unknowns] constructs an mfgSystem by building the \
structural equations (SignedFlows, Balance equations, HJ conditions, etc.) from \
the provided scenario and symbolic unknowns. makeSystem[s_scenario] automatically \
derives unknowns using makeUnknowns[s].";

SystemData::usage =
"SystemData[sys, key] returns the value associated with key in the system sys, or \
Missing[\"KeyAbsent\", key] if absent. SystemData[sys] returns the underlying \
Association.";

(* Structural Helpers *)

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

RoundValues::usage = "RoundValues[x] rounds numerical values in x to a standard precision (10^-10).";

Begin["`Private`"];

(* --- Structural Helpers Implementation --- *)

RoundValues[x_?NumberQ] :=
    Round[x, 10^-10]

RoundValues[Rule[a_, b_]] :=
    Rule[a, RoundValues[b]]

RoundValues[x_List] :=
    RoundValues /@ x

RoundValues[x_Association] :=
    RoundValues /@ x

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
        ];

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

(* --- Type predicate --- *)

mfgSystemQ[mfgSystem[_Association]] := True;
mfgSystemQ[_]                         := False;

(* --- Accessor --- *)

SystemData[mfgSystem[assoc_Association]]           := assoc;
SystemData[mfgSystem[assoc_Association], key_]     := Lookup[assoc, key, Missing["KeyAbsent", key]];

(* --- Constructors --- *)

makeSystem[s_?scenarioQ] := makeSystem[s, makeUnknowns[s]];

makeSystem[s_?scenarioQ, unk_?unknownsQ] :=
    Module[{model, topology, verticesList, adjacencyMatrix, entryVerticesFlows,
         exitVerticesCosts, graph, entryVertices, exitVertices, auxEntryVertices, 
         auxExitVertices, entryEdges, exitEdges, auxiliaryGraph, auxEdgeList, 
         edgeList, auxVerticesList, auxPairs, auxTriples, EntryDataAssociation,
         ExitCosts, js, us, jts, SignedFlows, SwitchingCosts, IneqJs, IneqJts,
         AltFlows, AltTransitionFlows, splittingPairs, EqBalanceSplittingFlows,
         BalanceSplittingFlows, NoDeadStarts, RuleBalanceGatheringFlows, BalanceGatheringFlows,
         EqBalanceGatheringFlows, EqEntryIn, RuleEntryValues, RuleEntryOut, RuleExitFlowsIn,
         RuleExitValues, EqValueAuxiliaryEdges, IneqSwitchingByVertex, AltOptCond,
         blockedIncomingPairs, blockedIncomingLookup, activeAuxTriples, auxTriplesByMiddle,
         scenarioAuxTriples, unknownAuxTriples,
         Nlhs, ModuleVars, ModuleVarsNames, Nrhs, costpluscurrents, EqGeneral,
         inAuxEntryPairs, outAuxEntryPairs, inAuxExitPairs, outAuxExitPairs, 
        pairs, halfPairs, consistentCosts, hamiltonian, alphaDefault, edgeAlpha, alphaAtEdge, edgeCost},
        
        model = ScenarioData[s, "Model"];
        topology = ScenarioData[s, "Topology"];
        hamiltonian = ScenarioData[s, "Hamiltonian"];
        If[!AssociationQ[hamiltonian],
            hamiltonian = <||>
        ];
        alphaDefault = Lookup[hamiltonian, "Alpha", 1];
        edgeAlpha = Lookup[hamiltonian, "EdgeAlpha", <||>];
        If[!AssociationQ[edgeAlpha],
            edgeAlpha = <||>
        ];
        alphaAtEdge[edge_List] :=
            Lookup[
                edgeAlpha,
                Key[edge],
                Lookup[edgeAlpha, Key[Reverse[edge]], alphaDefault]
            ];
        edgeCost[m_, edge_List] := m^alphaAtEdge[edge];
        
        If[!AssociationQ[topology],
            topology = BuildAuxiliaryTopology[model]
        ];
        
        If[topology === $Failed,
            Return[Failure["makeSystem", <|"Message" -> "Could not build topology from scenario."|>], Module]
        ];

        verticesList = model["Vertices List"];
        adjacencyMatrix = model["Adjacency Matrix"];
        entryVerticesFlows = model["Entrance Vertices and Flows"];
        exitVerticesCosts = model["Exit Vertices and Terminal Costs"];
        SwitchingCosts = model["Switching Costs"];

        graph = topology["Graph"];
        auxiliaryGraph = topology["AuxiliaryGraph"];
        auxEntryVertices = topology["AuxEntryVertices"];
        auxExitVertices = topology["AuxExitVertices"];
        entryEdges = topology["AuxEntryEdges"];
        exitEdges = topology["AuxExitEdges"];
        auxTriples = topology["AuxTriples"];
        
        halfPairs = topology["HalfPairs"];
        inAuxEntryPairs = topology["InAuxEntryPairs"];
        outAuxExitPairs = topology["OutAuxExitPairs"];
        inAuxExitPairs = topology["InAuxExitPairs"];
        outAuxEntryPairs = topology["OutAuxEntryPairs"];
        pairs = topology["Pairs"];
            
        auxEdgeList = EdgeList[auxiliaryGraph];
        edgeList = EdgeList[graph];
        auxVerticesList = VertexList[auxiliaryGraph];

        entryVertices = First /@ entryVerticesFlows;
        exitVertices = First /@ exitVerticesCosts;

        (* Prepare boundary data *)
        EntryDataAssociation = RoundValues @ AssociationThread[inAuxEntryPairs,
             Last /@ entryVerticesFlows];
        ExitCosts = AssociationThread[auxExitVertices, Last /@ exitVerticesCosts];

        (* Variables from unknown bundle *)
        js = UnknownsData[unk, "js"];
        jts = UnknownsData[unk, "jts"];
        us = UnknownsData[unk, "us"];
        auxPairs = UnknownsData[unk, "auxPairs"];
        unknownAuxTriples = UnknownsData[unk, "auxTriples"];

        (* Preserve canonical AuxTriples ordering from scenario topology:
           downstream transformations may depend on stable order. *)
        If[ListQ[unknownAuxTriples] && unknownAuxTriples =!= auxTriples,
            Return[
                Failure[
                    "makeSystem",
                    <|
                        "Message" -> "Mismatch between scenario topology AuxTriples and unknown bundle auxTriples; exact order must match.",
                        "ScenarioAuxTriplesLength" -> Length[auxTriples],
                        "UnknownsAuxTriplesLength" -> Length[unknownAuxTriples]
                    |>
                ],
                Module
            ]
        ];

        (* Signed flows *)
        SignedFlows = AssociationMap[j @@ # - j @@ Reverse @ #&, Join[
            inAuxEntryPairs, outAuxExitPairs, halfPairs]];

        consistentCosts = IsSwitchingCostConsistent[Normal @ SwitchingCosts];
        
        Which[
            consistentCosts === False,
                Message[mfgSystem::switchingcosts]
            ,
            consistentCosts =!= True,
                MFGPrint["Switching costs conditions are ", consistentCosts]
        ];

        IneqJs = And @@ (# >= 0& /@ js);
        IneqJts = And @@ (# >= 0& /@ jts);
        
        AltFlows = And @@ (AltFlowOp[j] /@ Join[
            inAuxEntryPairs, outAuxExitPairs, halfPairs]);
            
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
            
        splittingPairs = Join[inAuxEntryPairs, inAuxExitPairs, pairs];
        BalanceSplittingFlows = (j @@ # - Total[j @@@ FlowSplitting[auxTriples
            ][#]])& /@ splittingPairs;
        EqBalanceSplittingFlows = Simplify /@ (And @@ ((# == 0)& /@ BalanceSplittingFlows
            ));
            
        NoDeadStarts = Join[outAuxEntryPairs, outAuxExitPairs, pairs];
        RuleBalanceGatheringFlows = Association[(j @@ # -> Total[j @@@
             FlowGathering[auxTriples][#]])& /@ NoDeadStarts];
             
        BalanceGatheringFlows = ((-j @@ # + Total[j @@@ FlowGathering[
            auxTriples][#]])& /@ NoDeadStarts);
        EqBalanceGatheringFlows = Simplify /@ (And @@ (# == 0& /@ BalanceGatheringFlows
            ));
            
        EqEntryIn = (j @@ # == EntryDataAssociation[#])& /@ inAuxEntryPairs;
        
        RuleEntryOut = Association[(j @@ # -> 0)& /@ outAuxEntryPairs];
        RuleExitFlowsIn = Association[(j @@ # -> 0)& /@ inAuxExitPairs];
        
        blockedIncomingPairs = DeleteDuplicates @ Join[outAuxEntryPairs, inAuxExitPairs];
        blockedIncomingLookup =
            AssociationThread[blockedIncomingPairs, ConstantArray[True, Length[blockedIncomingPairs]]];
        activeAuxTriples = Select[
            auxTriples,
            !KeyExistsQ[blockedIncomingLookup, #[[{1, 2}]]] &
        ];
        auxTriplesByMiddle = GroupBy[auxTriples, #[[2]] &];
        
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
                    Lookup[auxTriplesByMiddle, #, {}],
                    !KeyExistsQ[blockedIncomingLookup, #[[{1, 2}]]] &
                ] & /@ verticesList;
        IneqSwitchingByVertex = And @@@ IneqSwitchingByVertex;
        
        AltOptCond = And @@ AltSwitch[j, u, SwitchingCosts] @@@ activeAuxTriples;
        
        Nlhs = Flatten[u @@ # - u @@ Reverse @ # + SignedFlows[#]& /@
             halfPairs];
        Nrhs = Flatten[SignedFlows[#] - Sign[SignedFlows[#]] edgeCost[SignedFlows[
            #], #]& /@ halfPairs];
            
        costpluscurrents = Table[Symbol["cpc" <> ToString[k]], {k, 1,
             Length @ edgeList}];
        EqGeneral = And @@ (MapThread[Equal, {Nlhs, costpluscurrents}
            ]);
        costpluscurrents = AssociationThread[costpluscurrents, Nrhs];
            
        ModuleVars = {graph, pairs, entryVertices, auxEntryVertices, 
            exitVertices,
            auxExitVertices, entryEdges, exitEdges, auxiliaryGraph,
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
             
        mfgSystem @ Join[
            model,
            AssociationThread[ModuleVarsNames, ModuleVars]
        ]
    ];

makeSystem[_] :=
    Failure["makeSystem", <|"Message" -> "makeSystem requires a scenario object."|>];

End[];
