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

BeginPackage["MFGraphs`"];

(* --- Public API declarations --- *)

mfgSystem::usage =
"mfgSystem[assoc] is the typed head for an MFG structural equation system. Use \
makeSystem to construct one; use SystemData to access keys.";

mfgBoundaryData::usage =
"mfgBoundaryData[assoc] is a typed record for boundary conditions and rules.";

mfgFlowData::usage =
"mfgFlowData[assoc] is a typed record for flow balance and non-negativity constraints.";

mfgComplementarityData::usage =
"mfgComplementarityData[assoc] is a typed record for complementarity conditions and switching costs.";

mfgHamiltonianData::usage =
"mfgHamiltonianData[assoc] is a typed record for Hamiltonian residuals and general equations.";

mfgSystem::switchingcosts = "Switching costs are inconsistent.";

mfgSystemQ::usage =
"mfgSystemQ[x] returns True if x is a typed mfgSystem[assoc_Association] object, \
False otherwise.";

mfgBoundaryDataQ::usage = "mfgBoundaryDataQ[x] returns True if x is a typed mfgBoundaryData object.";
mfgFlowDataQ::usage = "mfgFlowDataQ[x] returns True if x is a typed mfgFlowData object.";
mfgComplementarityDataQ::usage = "mfgComplementarityDataQ[x] returns True if x is a typed mfgComplementarityData object.";
mfgHamiltonianDataQ::usage = "mfgHamiltonianDataQ[x] returns True if x is a typed mfgHamiltonianData object.";

makeSystem::usage =
"makeSystem[s_scenario, unk_unknowns] constructs an mfgSystem by building the \
structural equations (SignedFlows, Balance equations, HJ conditions, etc.) from \
the provided scenario and symbolic unknowns. makeSystem[s_scenario] automatically \
derives unknowns using makeUnknowns[s].";

SystemData::usage =
"SystemData[sys, key] returns the value associated with key in the system sys, or \
Missing[\"KeyAbsent\", key] if absent. SystemData[sys] returns the underlying \
Association.";

SystemDataFlatten::usage =
"SystemDataFlatten[sys] returns a single flat Association containing all keys from \
all nested typed sub-records within the system. Useful for backward compatibility \
with legacy solvers.";

(* Structural Helpers *)

GetKirchhoffLinearSystem::usage = "GetKirchhoffLinearSystem[sys] returns the entry current vector, Kirchhoff matrix, and the variables in the order corresponding to the Kirchhoff matrix.";

GetKirchhoffMatrix::usage = "GetKirchhoffMatrix[sys] returns the entry current vector, Kirchhoff matrix, (critical congestion) cost function placeholder, and the variables in the order corresponding to the Kirchhoff matrix. The third slot is retained for backward compatibility.";

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

ConsistentSwitchingCosts[sc_Association][trip:{a_, b_, c_}] :=
    Module[{S = sc[trip], originTriples, conds},
        originTriples = DeleteCases[
            Select[Keys[sc], MatchQ[#, {a, b, _}] &],
            trip
        ];
        If[originTriples === {},
            True,
            conds = (S <= sc[#] + Lookup[sc, Key[{#[[3]], b, c}], Missing["KeyAbsent", {#[[3]], b, c}]]) & /@ originTriples;
            And @@ conds
        ]
    ];

ConsistentSwitchingCosts[sc_List][{a_, b_, c_} -> S_] :=
    ConsistentSwitchingCosts[Association[sc]][{a, b, c}];

IsSwitchingCostConsistent[switchingCosts_Association] :=
    Module[{triples, groupByAB, checkCost},
        triples = Keys[switchingCosts];
        groupByAB = GroupBy[triples, #[[;; 2]] &];
        checkCost[trip:{a_, b_, c_}] :=
            Module[{S = switchingCosts[trip], originTriples, conds},
                originTriples = DeleteCases[Lookup[groupByAB, Key[{a, b}], {}], trip];
                If[originTriples === {},
                    True,
                    conds = (S <= switchingCosts[#] + Lookup[switchingCosts, Key[{#[[3]], b, c}], Missing["KeyAbsent", {#[[3]], b, c}]]) & /@ originTriples;
                    And @@ conds
                ]
            ];
        And @@ Simplify[checkCost /@ triples]
    ];

IsSwitchingCostConsistent[switchingCosts_List] :=
    IsSwitchingCostConsistent[Association[switchingCosts]];

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

mfgBoundaryDataQ[mfgBoundaryData[_Association]] := True;
mfgBoundaryDataQ[_]                             := False;

mfgFlowDataQ[mfgFlowData[_Association]] := True;
mfgFlowDataQ[_]                         := False;

mfgComplementarityDataQ[mfgComplementarityData[_Association]] := True;
mfgComplementarityDataQ[_]                                     := False;

mfgHamiltonianDataQ[mfgHamiltonianData[_Association]] := True;
mfgHamiltonianDataQ[_]                                 := False;

(* --- Accessor --- *)

SystemData[mfgSystem[assoc_Association]]           := assoc;
SystemData[mfgSystem[assoc_Association], key_] :=
    Module[{val},
        val = Lookup[assoc, key, $Failed];
        If[val =!= $Failed, Return[val, Module]];
        
        (* Search in typed sub-records only *)
        Do[
            If[MatchQ[sub, (mfgBoundaryData|mfgFlowData|mfgComplementarityData|mfgHamiltonianData)[_Association]],
                val = Lookup[First[sub], key, $Failed];
                If[val =!= $Failed, Return[val, Module]]
            ],
            {sub, Values[assoc]}
        ];
        
        Missing["KeyAbsent", key]
    ];

SystemDataFlatten[mfgSystem[assoc_Association]] :=
    KeyDrop[
        Join[assoc, Sequence @@ (First /@ Select[Values[assoc],
            MatchQ[#, (mfgBoundaryData|mfgFlowData|mfgComplementarityData|mfgHamiltonianData)[_Association]]&])],
        {"BoundaryData", "FlowData", "ComplementarityData", "HamiltonianData"}
    ];

(* --- Modular Builders --- *)

BuildBoundaryData[s_?scenarioQ, topology_Association] :=
    Module[{model, entryVerticesFlows, exitVerticesCosts, inAuxEntryPairs,
         outAuxExitPairs, EntryDataAssociation,
         ExitCosts, EqEntryIn, RuleEntryIn, RuleEntryOut, RuleExitFlowsIn, RuleExitValues,
         auxExitVertices},
        model = ScenarioData[s, "Model"];
        entryVerticesFlows = model["Entrance Vertices and Flows"];
        exitVerticesCosts = model["Exit Vertices and Terminal Costs"];
        inAuxEntryPairs = topology["InAuxEntryPairs"];
        outAuxExitPairs = topology["OutAuxExitPairs"];
        auxExitVertices = topology["AuxExitVertices"];

        EntryDataAssociation = RoundValues @ AssociationThread[inAuxEntryPairs,
             Last /@ entryVerticesFlows];
        ExitCosts = AssociationThread[auxExitVertices, Last /@ exitVerticesCosts];
        
        EqEntryIn = (j @@ # == EntryDataAssociation[#])& /@ inAuxEntryPairs;
        RuleEntryIn = AssociationThread[j @@@ inAuxEntryPairs, Last /@ entryVerticesFlows];
        RuleEntryOut = <||>;
        RuleExitFlowsIn = <||>;
        
        RuleExitValues = AssociationThread[u @@@ outAuxExitPairs,
             Last /@ exitVerticesCosts];
        mfgBoundaryData @ <|
            "EntryDataAssociation" -> EntryDataAssociation,
            "ExitCosts" -> ExitCosts,
            "EqEntryIn" -> EqEntryIn,
            "RuleEntryIn" -> RuleEntryIn,
            "RuleEntryOut" -> RuleEntryOut,
            "RuleExitFlowsIn" -> RuleExitFlowsIn,
            "RuleExitValues" -> RuleExitValues
        |>
    ];

BuildFlowData[s_?scenarioQ, topology_Association, unk_?unknownsQ] :=
    Module[{inAuxEntryPairs, outAuxExitPairs, halfPairs, js, jts, pairs, auxTriples,
         SignedFlows, IneqJs, IneqJts, splittingPairs, BalanceSplittingFlows,
         EqBalanceSplittingFlows, NoDeadStarts, RuleBalanceGatheringFlows,
         BalanceGatheringFlows, EqBalanceGatheringFlows,
         splittingMaps, gatheringMaps, boundaryPairs},
        inAuxEntryPairs = topology["InAuxEntryPairs"];
        outAuxExitPairs = topology["OutAuxExitPairs"];
        halfPairs = topology["HalfPairs"];
        pairs = topology["Pairs"];
        auxTriples = topology["AuxTriples"];
        
        splittingMaps = GroupBy[auxTriples, #[[{1, 2}]] &];
        gatheringMaps = GroupBy[auxTriples, #[[{2, 3}]] &];

        js = UnknownsData[unk, "js"];
        jts = UnknownsData[unk, "jts"];
        boundaryPairs = Join[inAuxEntryPairs, outAuxExitPairs];

        SignedFlows = AssociationMap[
            If[MemberQ[boundaryPairs, #], j @@ #, j @@ # - j @@ Reverse @ #] &,
            Join[inAuxEntryPairs, outAuxExitPairs, halfPairs]
        ];

        IneqJs = And @@ (# >= 0& /@ js);
        IneqJts = And @@ (# >= 0& /@ jts);
        
        splittingPairs = Join[inAuxEntryPairs, pairs];
        BalanceSplittingFlows = (j @@ # - Total[j @@@ Lookup[splittingMaps, Key[#], {}]])& /@ splittingPairs;
        EqBalanceSplittingFlows = Simplify /@ (And @@ ((# == 0)& /@ BalanceSplittingFlows
            ));
            
        NoDeadStarts = Join[outAuxExitPairs, pairs];
        RuleBalanceGatheringFlows = Association[(j @@ # -> Total[j @@@
             Lookup[gatheringMaps, Key[#], {}]])& /@ NoDeadStarts];
             
        BalanceGatheringFlows = ((-j @@ # + Total[j @@@ Lookup[gatheringMaps, 
            Key[#], {}]])& /@ NoDeadStarts);
        EqBalanceGatheringFlows = Simplify /@ (And @@ (# == 0& /@ BalanceGatheringFlows
            ));

        mfgFlowData @ <|
            "SignedFlows" -> SignedFlows,
            "IneqJs" -> IneqJs,
            "IneqJts" -> IneqJts,
            "BalanceSplittingFlows" -> BalanceSplittingFlows,
            "EqBalanceSplittingFlows" -> EqBalanceSplittingFlows,
            "RuleBalanceGatheringFlows" -> RuleBalanceGatheringFlows,
            "BalanceGatheringFlows" -> BalanceGatheringFlows,
            "EqBalanceGatheringFlows" -> EqBalanceGatheringFlows,
            "splittingPairs" -> splittingPairs,
            "NoDeadStarts" -> NoDeadStarts
        |>
    ];

BuildComplementarityData[s_?scenarioQ, topology_Association, unk_?unknownsQ] :=
    Module[{model, verticesList, SwitchingCosts, inAuxEntryPairs, outAuxExitPairs,
         halfPairs, auxTriples, consistentCosts, AltFlows, AltTransitionFlows,
         activeAuxTriples, auxTriplesByMiddle, IneqSwitchingByVertex, AltOptCond,
         auxEntrySet, auxExitSet, interiorTripleQ},
        model = ScenarioData[s, "Model"];
        verticesList = model["Vertices List"];
        SwitchingCosts = model["Switching Costs"];
        inAuxEntryPairs = topology["InAuxEntryPairs"];
        outAuxExitPairs = topology["OutAuxExitPairs"];
        halfPairs = topology["HalfPairs"];
        auxTriples = topology["AuxTriples"];

        auxEntrySet = AssociationThread[topology["AuxEntryVertices"], True];
        auxExitSet  = AssociationThread[topology["AuxExitVertices"],  True];
        (* A triple is interior if neither endpoint is a boundary aux vertex.
           Reversed boundary triples are always 0, making AltFlowOp trivially true. *)
        interiorTripleQ[{r_, _, w_}] :=
            !KeyExistsQ[auxEntrySet, r] && !KeyExistsQ[auxExitSet, w] &&
            !KeyExistsQ[auxExitSet,  r] && !KeyExistsQ[auxEntrySet, w];

        consistentCosts = IsSwitchingCostConsistent[Normal @ SwitchingCosts];

        If[consistentCosts === False, Message[mfgSystem::switchingcosts]];

        AltFlows = And @@ (AltFlowOp[j] /@ halfPairs);
            
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
                    Select[auxTriples, interiorTripleQ[#] && OrderedQ[{First[#], Last[#]}]&])
            ];
            
        activeAuxTriples = auxTriples;
        auxTriplesByMiddle = GroupBy[auxTriples, #[[2]] &];
        
        IneqSwitchingByVertex =
            IneqSwitch[u, SwitchingCosts] @@@
                Lookup[auxTriplesByMiddle, #, {}] & /@ verticesList;
        IneqSwitchingByVertex = And @@@ IneqSwitchingByVertex;
        
        AltOptCond = And @@ AltSwitch[j, u, SwitchingCosts] @@@ activeAuxTriples;

        mfgComplementarityData @ <|
            "SwitchingCosts" -> SwitchingCosts,
            "consistentCosts" -> consistentCosts,
            "AltFlows" -> AltFlows,
            "AltTransitionFlows" -> AltTransitionFlows,
            "IneqSwitchingByVertex" -> IneqSwitchingByVertex,
            "AltOptCond" -> AltOptCond,
            "activeAuxTriples" -> activeAuxTriples
        |>
    ];

BuildHamiltonianData[s_?scenarioQ, topology_Association, flowData_mfgFlowData] :=
    Module[{hamiltonian, alphaDefault, edgeAlpha, alphaAtEdge, edgeCost, halfPairs,
         SignedFlows, Nlhs, Nrhs, costpluscurrents, EqGeneral},
        hamiltonian = ScenarioData[s, "Hamiltonian"];
        If[!AssociationQ[hamiltonian], hamiltonian = <||>];
        alphaDefault = Lookup[hamiltonian, "Alpha", 1];
        edgeAlpha = Lookup[hamiltonian, "EdgeAlpha", <||>];
        
        alphaAtEdge[edge_List] :=
            Lookup[
                edgeAlpha,
                Key[edge],
                Lookup[edgeAlpha, Key[Reverse[edge]], alphaDefault]
            ];
        edgeCost[m_, edge_List] := m^alphaAtEdge[edge];
        
        halfPairs = topology["HalfPairs"];
        SignedFlows = First[flowData]["SignedFlows"];

        Nlhs = Flatten[u @@ # - u @@ Reverse @ # + SignedFlows[#]& /@ halfPairs];
        Nrhs = Flatten[SignedFlows[#] - Sign[SignedFlows[#]] edgeCost[SignedFlows[#], #]& /@ halfPairs];

        costpluscurrents = Table[Symbol["MFGraphs`Private`cpc" <> ToString[k]], {k, 1, EdgeCount[topology["Graph"]]}];
        (* For alpha = 1 (critical congestion) the rhs is identically 0; use 0 directly. *)
        EqGeneral = And @@ MapThread[
            Equal[#1, If[alphaAtEdge[halfPairs[[#3]]] === 1, 0, #2]] &,
            {Nlhs, costpluscurrents, Range[Length[halfPairs]]}
        ];
        costpluscurrents = AssociationThread[costpluscurrents, Nrhs];

        mfgHamiltonianData @ <|
            "Nlhs" -> Nlhs,
            "Nrhs" -> Nrhs,
            "costpluscurrents" -> costpluscurrents,
            "EqGeneral" -> EqGeneral
        |>
    ];

(* --- Linear System Construction --- *)

GetKirchhoffLinearSystem[sys_] :=
    Module[{Kirchhoff, EqEntryIn, BalanceGatheringFlows, BalanceSplittingFlows,
         RuleEntryIn, RuleExitFlowsIn, RuleEntryOut, BM, KM, vars},
        EqEntryIn = SystemData[sys, "EqEntryIn"];
        BalanceGatheringFlows = SystemData[sys, "BalanceGatheringFlows"];
        BalanceSplittingFlows = SystemData[sys, "BalanceSplittingFlows"];
        RuleEntryIn = SystemData[sys, "RuleEntryIn"];
        RuleExitFlowsIn = SystemData[sys, "RuleExitFlowsIn"];
        RuleEntryOut = SystemData[sys, "RuleEntryOut"];

        (* Use defaults for missing data *)
        If[MissingQ[EqEntryIn], EqEntryIn = True];
        If[MissingQ[BalanceGatheringFlows], BalanceGatheringFlows = {}];
        If[MissingQ[BalanceSplittingFlows], BalanceSplittingFlows = {}];
        If[MissingQ[RuleEntryIn], RuleEntryIn = <||>];
        If[MissingQ[RuleExitFlowsIn], RuleExitFlowsIn = <||>];
        If[MissingQ[RuleEntryOut], RuleEntryOut = <||>];

        Kirchhoff = Join[
            If[Head[EqEntryIn] === List, EqEntryIn, {EqEntryIn}],
            (# == 0& /@ (BalanceGatheringFlows + BalanceSplittingFlows))
        ];
        Kirchhoff = Kirchhoff /. Join[Normal[RuleEntryIn], Normal[RuleExitFlowsIn], Normal[RuleEntryOut]];
        Kirchhoff = DeleteCases[Kirchhoff, True];
        Kirchhoff = Replace[Kirchhoff, False -> (0 == 1), {1}];
            
        vars = Select[Variables[Kirchhoff /. Equal -> List], MatchQ[#, j[_, _, _] | j[_, _]]&];
        If[Length[vars] === 0,
            {0, 0, {}}
            ,
            {BM, KM} = CoefficientArrays[Kirchhoff, vars];
            {-BM, KM, vars}
        ]
    ];

GetKirchhoffMatrix[sys_] :=
    Module[{B, KM, vars, cost, CCost},
        {B, KM, vars} = GetKirchhoffLinearSystem[sys];
        If[Length[vars] === 0,
            Return[{B, KM, <||>, vars}, Module]
        ];
        (* Legacy zero-cost placeholder *)
        cost = AssociationThread[vars, (0&) /@ vars];
        CCost = cost /@ vars /. MapThread[Rule, {vars, #}]&;
        {B, KM, CCost, vars}
    ];

(* --- Constructors --- *)

makeSystem[s_?scenarioQ] := makeSystem[s, makeUnknowns[s]];

makeSystem[s_?scenarioQ, unk_?unknownsQ] :=
    Module[{model, topology, graph, auxiliaryGraph, auxEdgeList, edgeList,
         auxVerticesList, auxTriples, unknownAuxTriples, boundaryData, flowData, compData, hamData,
         finalAssoc},
        
        model = ScenarioData[s, "Model"];
        topology = ScenarioData[s, "Topology"];
        If[!AssociationQ[topology], topology = BuildAuxiliaryTopology[model]];
        
        If[topology === $Failed,
            Return[Failure["makeSystem", <|"Message" -> "Could not build topology from scenario."|>], Module]
        ];

        (* Variables from scenario/unknowns needed for the top-level mfgSystem *)
        graph = topology["Graph"];
        auxiliaryGraph = topology["AuxiliaryGraph"];
        auxEdgeList = EdgeList[auxiliaryGraph];
        edgeList = EdgeList[graph];
        auxVerticesList = VertexList[auxiliaryGraph];
        auxTriples = topology["AuxTriples"];
        
        (* Preserve canonical AuxTriples ordering from scenario topology:
           downstream transformations may depend on stable order. *)
        unknownAuxTriples = UnknownsData[unk, "auxTriples"];
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

        (* Orchestrate Builders *)
        boundaryData = BuildBoundaryData[s, topology];
        flowData = BuildFlowData[s, topology, unk];
        compData = BuildComplementarityData[s, topology, unk];
        hamData = BuildHamiltonianData[s, topology, flowData];

        finalAssoc = Join[
            model,
            topology,
            <|
                "js" -> UnknownsData[unk, "js"],
                "us" -> UnknownsData[unk, "us"],
                "jts" -> UnknownsData[unk, "jts"],
                "edgeList" -> edgeList,
                "auxEdgeList" -> auxEdgeList,
                "auxVerticesList" -> auxVerticesList,
                "auxTriples" -> auxTriples,
                "BoundaryData" -> boundaryData,
                "FlowData" -> flowData,
                "ComplementarityData" -> compData,
                "HamiltonianData" -> hamData
            |>
        ];

        mfgSystem[finalAssoc]
    ];

makeSystem[_] :=
    Failure["makeSystem", <|"Message" -> "makeSystem requires a scenario object."|>];

End[];

EndPackage[];
