(* Wolfram Language package *)
(*
   systemTools.wl — typed system kernel for MFGraphs.

   Encapsulates the structural equations, inequalities, and substitution rules
   that represent a stationary Mean Field Game on a network.

   Lifecycle:
     makeSystem[s, unk]  →  builds system association + wraps in mfgSystem[...]
     makeSystem[s]       →  calls makeUnknowns[s] and delegates
     mfgSystemQ[x]       →  True iff x is a well-formed typed mfgSystem
     systemData[sys]     →  returns the underlying equation association
*)

BeginPackage["systemTools`", {"primitives`", "scenarioTools`", "unknownsTools`"}]

(* --- Public API declarations --- *)

mfgSystem::usage =
"mfgSystem[assoc] is the typed head for an MFG structural equation system. Use \
makeSystem to construct one; use systemData to access keys.";

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

systemData::usage =
"systemData[sys, key] returns the value associated with key in the system sys, or \
Missing[\"KeyAbsent\", key] if absent. systemData[sys] returns the underlying \
Association.";

systemDataFlatten::usage =
"systemDataFlatten[sys] returns a single flat Association containing all keys from \
all nested typed sub-records within the system. Useful for backward compatibility \
with legacy solvers.";

(* Structural Helpers *)

getKirchhoffLinearSystem::usage = "getKirchhoffLinearSystem[sys] returns the entry current vector, Kirchhoff matrix, and the variables in the order corresponding to the Kirchhoff matrix.";

getKirchhoffMatrix::usage = "getKirchhoffMatrix[sys] returns the entry current vector, Kirchhoff matrix, (critical congestion) cost function placeholder, and the variables in the order corresponding to the Kirchhoff matrix. The third slot is retained for backward compatibility.";

consistentSwitchingCosts::usage = "consistentSwitchingCosts[switchingcosts][{a,b,c}->S]
returns True if S, the cost of switching from the edge ab to cb, is smaller than any other combination,
such as, ab to bd and then from db to bc.
Returns the condition for this switching cost to satisfy the triangle inequality when S, and the other
switching costs too, does not have a numerical value.";

isSwitchingCostConsistent::usage = "isSwitchingCostConsistent[List of switching costs] is True if all switching costs satisfy the triangle inequality. If some switching costs are symbolic, then it returns the consistency conditions."

altFlowOp::usage = "altFlowOp[j][list] returns the alternative: j@@list ==0 || j@@Reverse@list ==0.";

flowSplitting::usage = "flowSplitting[AT][UndirectedEdge[a, b]] returns the splitting that start with {a,b}.";

flowGathering::usage = "flowGathering[auxTriples_List][x_] returns the triples that end with x.";

ineqSwitch::usage = "ineqSwitch[u, switchingCosts][v, e1, e2] returns the optimality condition at the vertex v related to switching from e1 to e2. Namely,
u[v, e1] <= u[e2, e1] + switchingCosts[{v, e1, e2}]"

altSwitch::usage = "altSwitch[j, u, switchingCosts][v, e1, e2] returns the complementarity condition:
(j[v, e1, e2] == 0) || (u[v, e1] == u[e2, e1] + switchingCosts[{v, e1, e2}])"

roundValues::usage = "roundValues[x] rounds numerical values in x to a standard precision (10^-10).";

Begin["`Private`"]

(* --- Structural Helpers Implementation --- *)

roundValues[x_?NumberQ] :=
    Round[x, 10^-10]

roundValues[Rule[a_, b_]] :=
    Rule[a, roundValues[b]]

roundValues[x_List] :=
    roundValues /@ x

roundValues[x_Association] :=
    roundValues /@ x

exactBoundaryValue[val_] /; ExactNumberQ[val] := val;

exactBoundaryValue[val_?NumericQ] :=
    Module[{nval = N[val]},
        If[TrueQ[Im[nval] == 0],
            Rationalize[val, 0],
            Failure["buildBoundaryData", <|
                "Message" -> "Boundary values must be real numeric values.",
                "Value" -> val
            |>]
        ]
    ];

exactBoundaryValue[val_] :=
    Failure["buildBoundaryData", <|
        "Message" -> "Boundary values must be numeric values.",
        "Value" -> val
    |>];

exactBoundaryValues[vals_List] :=
    Module[{exactVals = exactBoundaryValue /@ vals, firstFailure},
        firstFailure = FirstCase[exactVals, _Failure, Missing["NotFound"]];
        If[FailureQ[firstFailure], firstFailure, exactVals]
    ];

consistentSwitchingCosts[sc_Association][trip:{a_, b_, c_}] :=
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

consistentSwitchingCosts[sc_List][{a_, b_, c_} -> S_] :=
    consistentSwitchingCosts[Association[sc]][{a, b, c}];

isSwitchingCostConsistent[switchingCosts_Association] :=
    Module[{triples, groupByAB, checkCost},
        triples = Keys[switchingCosts];
        groupByAB = GroupBy[triples, #[[;; 2]] &];
        checkCost[trip:{a_, b_, c_}] :=
            Module[{S = switchingCosts[trip], originTriples, conds},
                originTriples = DeleteCases[Lookup[groupByAB, Key[{a, b}], {}], trip];
                If[originTriples === {},
                    True,
                    conds = Table[
                        With[{leg2 = Lookup[switchingCosts, Key[{trip2[[3]], b, c}], Missing["Impossible"]]},
                            If[MissingQ[leg2], True, S <= switchingCosts[trip2] + leg2]
                        ],
                        {trip2, originTriples}
                    ];
                    And @@ conds
                ]
            ];
        And @@ Simplify[checkCost /@ triples]
    ];

isSwitchingCostConsistent[switchingCosts_List] :=
    isSwitchingCostConsistent[Association[switchingCosts]];

altFlowOp[j_][list_] :=
    j @@ list == 0 || j @@ Reverse @ list == 0;

flowSplitting[auxTriples_List][x_] :=
    Select[auxTriples, MatchQ[#, {Sequence @@ x, __}]&];

flowGathering[auxTriples_List][x_] :=
    Select[auxTriples, MatchQ[#, {__, Sequence @@ x}]&]

ineqSwitch[u_, switching_Association][r_, i_, w_] :=
    u[r, i] <= u[w, i] + switching[{r, i, w}];

altSwitch[j_, u_, switching_][r_, i_, w_] :=
    (j[r, i, w] == 0) || (u[r, i] == u[w, i] + switching[{r, i, w}]);

(* --- Type predicates --- *)

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

systemData[mfgSystem[assoc_Association]]           := assoc;
systemData[mfgSystem[assoc_Association], key_] :=
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

systemDataFlatten[mfgSystem[assoc_Association]] :=
    KeyDrop[
        Join[assoc, Sequence @@ (First /@ Select[Values[assoc],
            MatchQ[#, (mfgBoundaryData|mfgFlowData|mfgComplementarityData|mfgHamiltonianData)[_Association]]&])],
        {"BoundaryData", "FlowData", "ComplementarityData", "HamiltonianData"}
    ];

(* --- Modular Builders --- *)

buildBoundaryData[s_?scenarioQ, topology_Association] :=
    Module[{model, entryVerticesFlows, exitVerticesCosts, inAuxEntryPairs,
         outAuxExitPairs, entryDataAssoc,
         exitCosts, eqEntryIn, ruleEntryIn, ruleEntryOut, ruleExitFlowsIn, ruleExitValues,
         auxExitVertices, entryValues, exitValues},
        model = scenarioData[s, "Model"];
        entryVerticesFlows = model["Entries"];
        exitVerticesCosts  = model["Exits"];
        inAuxEntryPairs = topology["InAuxEntryPairs"];
        outAuxExitPairs = topology["OutAuxExitPairs"];
        auxExitVertices = topology["AuxExitVertices"];

        entryValues = exactBoundaryValues[Last /@ entryVerticesFlows];
        If[FailureQ[entryValues], Return[entryValues, Module]];
        exitValues = exactBoundaryValues[Last /@ exitVerticesCosts];
        If[FailureQ[exitValues], Return[exitValues, Module]];

        entryDataAssoc = AssociationThread[inAuxEntryPairs, entryValues];
        exitCosts      = AssociationThread[auxExitVertices, exitValues];

        eqEntryIn      = (j @@ # == entryDataAssoc[#])& /@ inAuxEntryPairs;
        ruleEntryIn    = AssociationThread[j @@@ inAuxEntryPairs, entryValues];
        ruleEntryOut   = <||>;
        ruleExitFlowsIn = <||>;

        ruleExitValues = AssociationThread[u @@@ (Reverse /@ outAuxExitPairs), exitValues];

        mfgBoundaryData @ <|
            "EntryDataAssociation" -> entryDataAssoc,
            "ExitCosts"            -> exitCosts,
            "EqEntryIn"            -> eqEntryIn,
            "RuleEntryIn"          -> ruleEntryIn,
            "RuleEntryOut"         -> ruleEntryOut,
            "RuleExitFlowsIn"      -> ruleExitFlowsIn,
            "RuleExitValues"       -> ruleExitValues
        |>
    ];

buildFlowData[s_?scenarioQ, topology_Association, unk_?unknownsQ] :=
    Module[{inAuxEntryPairs, outAuxExitPairs, halfPairs, js, jts, pairs, auxTriples,
         signedFlows, ineqJs, ineqJts, splittingPairs, balanceSplittingFlows,
         eqBalanceSplittingFlows, noDeadStarts, ruleBalanceGatheringFlows,
         balanceGatheringFlows, eqBalanceGatheringFlows,
         splittingMaps, gatheringMaps, boundaryPairs},
        inAuxEntryPairs = topology["InAuxEntryPairs"];
        outAuxExitPairs = topology["OutAuxExitPairs"];
        halfPairs   = topology["HalfPairs"];
        pairs       = topology["Pairs"];
        auxTriples  = topology["AuxTriples"];

        splittingMaps = GroupBy[auxTriples, #[[{1, 2}]] &];
        gatheringMaps = GroupBy[auxTriples, #[[{2, 3}]] &];

        js = unknownsData[unk, "Js"];
        jts = unknownsData[unk, "Jts"];
        boundaryPairs = Join[inAuxEntryPairs, outAuxExitPairs];

        signedFlows = AssociationMap[
            If[MemberQ[boundaryPairs, #], j @@ #, j @@ # - j @@ Reverse @ #] &,
            Join[inAuxEntryPairs, outAuxExitPairs, halfPairs]
        ];

        ineqJs  = And @@ (# >= 0& /@ js);
        ineqJts = And @@ (# >= 0& /@ jts);

        splittingPairs = Join[inAuxEntryPairs, pairs];
        balanceSplittingFlows = (j @@ # - Total[j @@@ Lookup[splittingMaps, Key[#], {}]])& /@ splittingPairs;
        eqBalanceSplittingFlows = Simplify /@ (And @@ ((# == 0)& /@ balanceSplittingFlows));

        noDeadStarts = Join[outAuxExitPairs, pairs];
        ruleBalanceGatheringFlows = Association[(j @@ # -> Total[j @@@
             Lookup[gatheringMaps, Key[#], {}]])& /@ noDeadStarts];

        balanceGatheringFlows = ((-j @@ # + Total[j @@@ Lookup[gatheringMaps,
            Key[#], {}]])& /@ noDeadStarts);
        eqBalanceGatheringFlows = Simplify /@ (And @@ (# == 0& /@ balanceGatheringFlows));

        mfgFlowData @ <|
            "SignedFlows"              -> signedFlows,
            "IneqJs"                   -> ineqJs,
            "IneqJts"                  -> ineqJts,
            "BalanceSplittingFlows"    -> balanceSplittingFlows,
            "EqBalanceSplittingFlows"  -> eqBalanceSplittingFlows,
            "RuleBalanceGatheringFlows" -> ruleBalanceGatheringFlows,
            "BalanceGatheringFlows"    -> balanceGatheringFlows,
            "EqBalanceGatheringFlows"  -> eqBalanceGatheringFlows,
            "SplittingPairs"           -> splittingPairs,
            "NoDeadStarts"             -> noDeadStarts
        |>
    ];

buildComplementarityData[s_?scenarioQ, topology_Association, unk_?unknownsQ] :=
    Module[{model, verticesList, switchingCosts, inAuxEntryPairs, outAuxExitPairs,
         halfPairs, auxTriples, consistCosts, altFlows, altTransitionFlows,
         activeTriples, auxTriplesByMiddle, ineqSwitchingByVertex, altOptCond,
         auxEntrySet, auxExitSet, interiorTripleQ},
        model = scenarioData[s, "Model"];
        verticesList  = model["Vertices"];
        switchingCosts = model["Switching"];
        inAuxEntryPairs = topology["InAuxEntryPairs"];
        outAuxExitPairs = topology["OutAuxExitPairs"];
        halfPairs  = topology["HalfPairs"];
        auxTriples = topology["AuxTriples"];

        auxEntrySet = AssociationThread[topology["AuxEntryVertices"], True];
        auxExitSet  = AssociationThread[topology["AuxExitVertices"],  True];
        interiorTripleQ[{r_, _, w_}] :=
            !KeyExistsQ[auxEntrySet, r] && !KeyExistsQ[auxExitSet, w] &&
            !KeyExistsQ[auxExitSet,  r] && !KeyExistsQ[auxEntrySet, w];

        consistCosts = isSwitchingCostConsistent[Normal @ switchingCosts];

        If[consistCosts === False, Message[mfgSystem::switchingcosts]];

        altFlows = And @@ (altFlowOp[j] /@ halfPairs);

        altTransitionFlows =
            If[consistCosts === False,
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
                And @@ (altFlowOp[j] /@
                    Select[auxTriples, interiorTripleQ[#] && OrderedQ[{First[#], Last[#]}]&])
            ];

        activeTriples = auxTriples;
        auxTriplesByMiddle = GroupBy[auxTriples, #[[2]] &];

        ineqSwitchingByVertex =
            ineqSwitch[u, switchingCosts] @@@
                Lookup[auxTriplesByMiddle, #, {}] & /@ verticesList;
        ineqSwitchingByVertex = And @@ Flatten[ineqSwitchingByVertex];

        altOptCond = And @@ altSwitch[j, u, switchingCosts] @@@ activeTriples;

        mfgComplementarityData @ <|
            "SwitchingCosts"         -> switchingCosts,
            "ConsistentCosts"        -> consistCosts,
            "AltFlows"               -> altFlows,
            "AltTransitionFlows"     -> altTransitionFlows,
            "IneqSwitchingByVertex"  -> ineqSwitchingByVertex,
            "AltOptCond"             -> altOptCond,
            "ActiveTriples"          -> activeTriples
        |>
    ];

buildHamiltonianData[s_?scenarioQ, topology_Association, flowData_mfgFlowData] :=
    Module[{hamiltonian, alphaDefault, edgeAlpha, alphaAtEdge, edgeCost, halfPairs,
         signedFlows, nlhs, nrhs, costCurrents, eqGeneral},
        hamiltonian = scenarioData[s, "Hamiltonian"];
        If[!AssociationQ[hamiltonian], hamiltonian = <||>];
        alphaDefault = Lookup[hamiltonian, "Alpha", 1];
        edgeAlpha    = Lookup[hamiltonian, "EdgeAlpha", <||>];

        alphaAtEdge[edge_List] :=
            Lookup[
                edgeAlpha,
                Key[edge],
                Lookup[edgeAlpha, Key[Reverse[edge]], alphaDefault]
            ];
        edgeCost[m_, edge_List] := m^alphaAtEdge[edge];

        halfPairs   = topology["HalfPairs"];
        signedFlows = First[flowData]["SignedFlows"];

        nlhs = Flatten[u @@ # - u @@ Reverse @ # + signedFlows[#]& /@ halfPairs];
        nrhs = Flatten[signedFlows[#] - Sign[signedFlows[#]] edgeCost[signedFlows[#], #]& /@ halfPairs];

        costCurrents = Table[Symbol["systemTools`Private`cpc" <> ToString[k]], {k, 1, EdgeCount[topology["Graph"]]}];
        eqGeneral = And @@ MapThread[
            Equal[#1, If[alphaAtEdge[halfPairs[[#3]]] === 1, 0, #2]] &,
            {nlhs, costCurrents, Range[Length[halfPairs]]}
        ];
        costCurrents = AssociationThread[costCurrents, nrhs];

        mfgHamiltonianData @ <|
            "Nlhs"         -> nlhs,
            "Nrhs"         -> nrhs,
            "CostCurrents" -> costCurrents,
            "EqGeneral"    -> eqGeneral
        |>
    ];

(* --- Linear System Construction --- *)

getKirchhoffLinearSystem[sys_] :=
    Module[{kirchhoff, eqEntryIn, balanceGatheringFlows, balanceSplittingFlows,
         ruleEntryIn, ruleExitFlowsIn, ruleEntryOut, bm, km, vars},
        eqEntryIn            = systemData[sys, "EqEntryIn"];
        balanceGatheringFlows = systemData[sys, "BalanceGatheringFlows"];
        balanceSplittingFlows = systemData[sys, "BalanceSplittingFlows"];
        ruleEntryIn          = systemData[sys, "RuleEntryIn"];
        ruleExitFlowsIn      = systemData[sys, "RuleExitFlowsIn"];
        ruleEntryOut         = systemData[sys, "RuleEntryOut"];

        If[MissingQ[eqEntryIn],            eqEntryIn = True];
        If[MissingQ[balanceGatheringFlows], balanceGatheringFlows = {}];
        If[MissingQ[balanceSplittingFlows], balanceSplittingFlows = {}];
        If[MissingQ[ruleEntryIn],          ruleEntryIn = <||>];
        If[MissingQ[ruleExitFlowsIn],      ruleExitFlowsIn = <||>];
        If[MissingQ[ruleEntryOut],         ruleEntryOut = <||>];

        kirchhoff = Join[
            If[Head[eqEntryIn] === List, eqEntryIn, {eqEntryIn}],
            (# == 0& /@ (balanceGatheringFlows + balanceSplittingFlows))
        ];
        kirchhoff = kirchhoff /. Join[Normal[ruleEntryIn], Normal[ruleExitFlowsIn], Normal[ruleEntryOut]];
        kirchhoff = DeleteCases[kirchhoff, True];
        kirchhoff = Replace[kirchhoff, False -> (0 == 1), {1}];

        vars = Select[Variables[kirchhoff /. Equal -> List], MatchQ[#, j[_, _, _] | j[_, _]]&];
        If[Length[vars] === 0,
            {0, 0, {}}
            ,
            {bm, km} = CoefficientArrays[kirchhoff, vars];
            {-bm, km, vars}
        ]
    ];

getKirchhoffMatrix[sys_] :=
    Module[{b, km, vars, cost, cCost},
        {b, km, vars} = getKirchhoffLinearSystem[sys];
        If[Length[vars] === 0,
            Return[{b, km, <||>, vars}, Module]
        ];
        cost  = AssociationThread[vars, (0&) /@ vars];
        cCost = cost /@ vars /. MapThread[Rule, {vars, #}]&;
        {b, km, cCost, vars}
    ];

(* --- Constructors --- *)

makeSystem[s_?scenarioQ] := makeSystem[s, makeUnknowns[s]];

makeSystem[s_?scenarioQ, unk_?unknownsQ] :=
    Module[{model, topology, graph, auxiliaryGraph, auxEdges, edges,
         auxVertices, auxTriples, unknownAuxTriples, boundaryData, flowData, compData, hamData,
         finalAssoc},

        model    = scenarioData[s, "Model"];
        topology = scenarioData[s, "Topology"];
        If[!AssociationQ[topology], topology = buildAuxiliaryTopology[model]];

        If[topology === $Failed,
            Return[Failure["makeSystem", <|"Message" -> "Could not build topology from scenario."|>], Module]
        ];

        graph          = topology["Graph"];
        auxiliaryGraph = topology["AuxiliaryGraph"];
        auxEdges       = EdgeList[auxiliaryGraph];
        edges          = EdgeList[graph];
        auxVertices    = VertexList[auxiliaryGraph];
        auxTriples     = topology["AuxTriples"];

        (* Preserve canonical AuxTriples ordering from scenario topology *)
        unknownAuxTriples = unknownsData[unk, "AuxTriples"];
        If[ListQ[unknownAuxTriples] && unknownAuxTriples =!= auxTriples,
            Return[
                Failure[
                    "makeSystem",
                    <|
                        "Message" -> "Mismatch between scenario topology AuxTriples and unknown bundle AuxTriples; exact order must match.",
                        "ScenarioAuxTriplesLength" -> Length[auxTriples],
                        "UnknownsAuxTriplesLength" -> Length[unknownAuxTriples]
                    |>
                ],
                Module
            ]
        ];

        boundaryData = buildBoundaryData[s, topology];
        If[FailureQ[boundaryData],
            Return[
                Failure["makeSystem", <|
                    "Message" -> Lookup[boundaryData[[2]], "Message", "Failed to build boundary data."],
                    "Cause" -> boundaryData
                |>],
                Module
            ]
        ];
        flowData = buildFlowData[s, topology, unk];
        compData = buildComplementarityData[s, topology, unk];
        hamData  = buildHamiltonianData[s, topology, flowData];

        finalAssoc = Join[
            model,
            topology,
            <|
                "Js"          -> unknownsData[unk, "Js"],
                "Us"          -> unknownsData[unk, "Us"],
                "Jts"         -> unknownsData[unk, "Jts"],
                "Edges"       -> edges,
                "AuxEdges"    -> auxEdges,
                "AuxVertices" -> auxVertices,
                "AuxTriples"  -> auxTriples,
                "BoundaryData"       -> boundaryData,
                "FlowData"           -> flowData,
                "ComplementarityData" -> compData,
                "HamiltonianData"    -> hamData
            |>
        ];

        mfgSystem[finalAssoc]
    ];

makeSystem[_] :=
    Failure["makeSystem", <|"Message" -> "makeSystem requires a scenario object."|>];

End[]

EndPackage[]
