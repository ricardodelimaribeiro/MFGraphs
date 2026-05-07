(* Wolfram Language package *)
(*
   systemTools.wl — typed system kernel for MFGraphs.

   Encapsulates the structural equations, inequalities, and substitution rules
   that represent a stationary Mean Field Game on a network.

   Lifecycle:
     makeSystem[s, unk]  →  builds system association + wraps in mfgSystem[...]
     makeSystem[s]       →  calls makeSymbolicUnknowns[s] and delegates
     mfgSystemQ[x]       →  True iff x is a well-formed typed mfgSystem
     systemData[sys]     →  returns the underlying equation association
*)

BeginPackage["systemTools`", {"primitives`", "utilities`", "scenarioTools`", "unknownsTools`"}]

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
"makeSystem[s_scenario, unk_symbolicUnknowns] constructs an mfgSystem by building the \
structural equations (SignedFlows, Balance equations, HJ conditions, etc.) from \
the provided scenario and exact symbolic unknown bundle. makeSystem[s_scenario] \
automatically derives symbolicUnknowns using makeSymbolicUnknowns[s].";

buildBoundaryData::usage =
"buildBoundaryData[s, topology] builds typed boundary equations, entry rules, exit inequalities, and \
entry/exit metadata. IneqExitValues stores upper-bound inequalities u[auxExit,N]<=cost for each exit; \
AltExitCond stores the complementarity conditions j[N,auxExit]==0||u[auxExit,N]==cost.";

buildFlowData::usage =
"buildFlowData[s, topology, unk] builds typed flow-balance equations and non-negativity constraints.";

buildComplementarityData::usage =
"buildComplementarityData[s, topology, unk] builds typed complementarity alternatives and switching inequalities.";

buildHamiltonianData::usage =
"buildHamiltonianData[s, topology, flowData] builds typed Hamiltonian residual equations for the system. \
EqGeneral encodes the edge-level HJB equation u[a,b]-u[b,a]+j[a,b]-j[b,a] = nrhs unconditionally for every \
undirected edge {a,b}, where nrhs = 0 for Alpha==1 and m - Sign[m] m^alpha otherwise. \
The equation is enforced even when net flow is zero (zero-flow edges force u[a,b]=u[b,a]), which propagates \
through switching inequalities to pin value variables at bypassed exit nodes to their terminal cost. \
Current system construction uses Alpha/EdgeAlpha; V/G/EdgeV/EdgeG are preserved on scenarios for future work.";

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

altFlowOp::usage = "altFlowOp[j][list] returns the alternative: j@@list ==0 || j@@Reverse@list ==0.";

flowSplitting::usage = "flowSplitting[AT][UndirectedEdge[a, b]] returns the splitting that start with {a,b}.";

flowGathering::usage = "flowGathering[auxTriples_List][x_] returns the triples that end with x.";

switchingCostLookup::usage =
"switchingCostLookup[sc] returns a function f such that f[r, i, w] gives the switching cost \
for the transition from e_{r,i} to e_{i,w}, defaulting to 0 if absent.";

ineqSwitch::usage = "ineqSwitch[u, switchingCosts][r, i, w] returns the optimality condition at junction i for the transition from r to w. Namely,
u[w, i] + switchingCosts[r, i, w] - u[r, i] >= 0.";

altSwitch::usage = "altSwitch[j, u, switchingCosts][r, i, w] returns the complementarity condition at junction i for the transition from r to w:
(j[r, i, w] == 0) || (u[w, i] + switchingCosts[r, i, w] - u[r, i] == 0).";

(* The following function is declared and implemented in utilities.wl:
   - roundValues
*)

Begin["`Private`"]

(* The following functions are declared and implemented in utilities.wl:
   - exactBoundaryValue
   - exactBoundaryValues
*)

interiorTripleQ::usage =
"interiorTripleQ[triple, auxEntrySet, auxExitSet] returns True if the triple \
does not involve any auxiliary entry or exit vertices.";

(* --- Structural Helpers Implementation --- *)

altFlowOp[j_][list_] :=
    j @@ list == 0 || j @@ Reverse @ list == 0;

flowSplitting[auxTriples_List][x_] :=
    Select[auxTriples, MatchQ[#, {Sequence @@ x, __}]&];

flowGathering[auxTriples_List][x_] :=
    Select[auxTriples, MatchQ[#, {__, Sequence @@ x}]&]

switchingCostLookup[sc_Association][r_, i_, w_] := Lookup[sc, Key[{r, i, w}], 0]

ineqSwitch[u_, switching_][r_, i_, w_] :=
    Module[{cost = switching[r, i, w]},
        u[w, i] + If[MatchQ[cost, _Function], cost[j[r, i, w]], cost] - u[r, i] >= 0
    ];

altSwitch[j_, u_, switching_][r_, i_, w_] :=
    Module[{cost = switching[r, i, w]},
        (j[r, i, w] == 0) || (u[w, i] + If[MatchQ[cost, _Function], cost[j[r, i, w]], cost] - u[r, i] == 0)
    ];

interiorTripleQ[triple:{r_, _, w_}, auxEntrySet_, auxExitSet_] :=
    !KeyExistsQ[auxEntrySet, r] && !KeyExistsQ[auxExitSet, w] &&
    !KeyExistsQ[auxExitSet,  r] && !KeyExistsQ[auxEntrySet, w];

(* --- Type predicates --- *)

mfgSystemQ[x_] := mfgTypedQ[x, mfgSystem];

mfgBoundaryDataQ[x_] := mfgTypedQ[x, mfgBoundaryData];

mfgFlowDataQ[x_] := mfgTypedQ[x, mfgFlowData];

mfgComplementarityDataQ[x_] := mfgTypedQ[x, mfgComplementarityData];

mfgHamiltonianDataQ[x_] := mfgTypedQ[x, mfgHamiltonianData];

(* --- Accessor --- *)

systemData[sys_] := mfgData[sys];
systemData[sys_, key_] :=
    Module[{val, assoc},
        assoc = mfgData[sys];
        val = Lookup[assoc, key, $Failed];
        If[val =!= $Failed, Return[val, Module]];

        (* Search in typed sub-records only *)
        Do[
            If[mfgTypedQ[sub, mfgBoundaryData] || mfgTypedQ[sub, mfgFlowData] || 
               mfgTypedQ[sub, mfgComplementarityData] || mfgTypedQ[sub, mfgHamiltonianData],
                val = Lookup[mfgData[sub], key, $Failed];
                If[val =!= $Failed, Return[val, Module]]
            ],
            {sub, Values[assoc]}
        ];

        Missing["KeyAbsent", key]
    ];

systemDataFlatten[sys_] :=
    Module[{assoc = mfgData[sys]},
        KeyDrop[
            Join[assoc, Sequence @@ (mfgData /@ Select[Values[assoc],
                mfgTypedQ[#, mfgBoundaryData] || mfgTypedQ[#, mfgFlowData] || 
                mfgTypedQ[#, mfgComplementarityData] || mfgTypedQ[#, mfgHamiltonianData] &])],
            {"BoundaryData", "FlowData", "ComplementarityData", "HamiltonianData"}
        ]
    ];

(* --- Modular Builders --- *)

buildBoundaryData[s_?scenarioQ, topology_Association] :=
    Module[{model, entryVerticesFlows, exitVerticesCosts, inAuxEntryPairs,
         outAuxExitPairs, entryDataAssoc,
         exitCosts, eqEntryIn, ruleEntryIn, ruleEntryOut,
         ineqExitValues, altExitCond,
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

        eqEntryIn   = (j @@ # == entryDataAssoc[#])& /@ inAuxEntryPairs;
        ruleEntryIn = AssociationThread[j @@@ inAuxEntryPairs, entryValues];
        ruleEntryOut = <||>;

        ineqExitValues = MapThread[(u @@ Reverse[#1]) <= #2 &, {outAuxExitPairs, exitValues}];
        altExitCond    = MapThread[(j @@ #1 == 0) || ((u @@ Reverse[#1]) == #2) &,
                             {outAuxExitPairs, exitValues}];

        mfgBoundaryData @ <|
            "EntryDataAssociation" -> entryDataAssoc,
            "ExitCosts"            -> exitCosts,
            "EqEntryIn"            -> eqEntryIn,
            "RuleEntryIn"          -> ruleEntryIn,
            "RuleEntryOut"         -> ruleEntryOut,
            "IneqExitValues"       -> ineqExitValues,
            "AltExitCond"          -> altExitCond
        |>
    ];

buildFlowData[s_?scenarioQ, topology_Association, unk_?symbolicUnknownsQ] :=
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

        js = symbolicUnknownsData[unk, "Js"];
        jts = symbolicUnknownsData[unk, "Jts"];
        boundaryPairs = Join[inAuxEntryPairs, outAuxExitPairs];

        signedFlows = AssociationMap[
            If[MemberQ[boundaryPairs, #], j @@ #, j @@ # - j @@ Reverse @ #] &,
            Join[inAuxEntryPairs, outAuxExitPairs, halfPairs]
        ];

        ineqJs  = # >= 0& /@ js;
        ineqJts = # >= 0& /@ jts;

        splittingPairs = Join[inAuxEntryPairs, pairs];
        balanceSplittingFlows = (j @@ # - Total[j @@@ Lookup[splittingMaps, Key[#], {}]])& /@ splittingPairs;
        eqBalanceSplittingFlows = Simplify /@ (And @@ ((# == 0)& /@ balanceSplittingFlows));

        noDeadStarts = Join[pairs, outAuxExitPairs];
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

buildComplementarityData[s_?scenarioQ, topology_Association, unk_?symbolicUnknownsQ] :=
    Module[{model, verticesList, switchingCosts, inAuxEntryPairs, outAuxExitPairs,
         halfPairs, auxTriples, consistCosts, altFlows, altTransitionFlows,
         activeTriples, auxTriplesByMiddle, ineqSwitchingByVertex, altOptCond,
         auxEntrySet, auxExitSet,
         zeroSwitchUEqualities},
        model = scenarioData[s, "Model"];
        verticesList  = model["Vertices"];
        switchingCosts = model["Switching"];
        inAuxEntryPairs = topology["InAuxEntryPairs"];
        outAuxExitPairs = topology["OutAuxExitPairs"];
        halfPairs  = topology["HalfPairs"];
        auxTriples = topology["AuxTriples"];

        auxEntrySet = AssociationThread[topology["AuxEntryVertices"], True];
        auxExitSet  = AssociationThread[topology["AuxExitVertices"],  True];

        (* Switching Costs are enforced consistent during makeScenario *)
        consistCosts = True;

        altFlows = altFlowOp[j] /@ halfPairs;

        altTransitionFlows =
            altFlowOp[j] /@
                Select[auxTriples, interiorTripleQ[#, auxEntrySet, auxExitSet] && OrderedQ[{First[#], Last[#]}]&];

        activeTriples = auxTriples;
        auxTriplesByMiddle = GroupBy[auxTriples, #[[2]] &];

        With[{scFn = switchingCostLookup[switchingCosts]},
            ineqSwitchingByVertex =
                ineqSwitch[u, scFn] @@@
                    Lookup[auxTriplesByMiddle, #, {}] & /@ verticesList;
            ineqSwitchingByVertex = Flatten[ineqSwitchingByVertex];

            (* Derive u-equality rules forced by zero-cost symmetric triple pairs.
               Theorem: {r,v,w} and {w,v,r} both with cost 0 produce
               u[w,v]>=u[r,v] and u[r,v]>=u[w,v] from IneqSwitching; antisymmetry
               of >= forces u[r,v]=u[w,v]. Only real (non-aux) vertex pairs qualify. *)
            zeroSwitchUEqualities = Flatten @ Map[
                Function[v,
                    Module[{triples, realZeroTriples, tripleSet, pairs, edges, components},
                        triples = Lookup[auxTriplesByMiddle, v, {}];
                        realZeroTriples = Select[triples, interiorTripleQ[#, auxEntrySet, auxExitSet] && scFn @@ # == 0 &];
                        If[realZeroTriples === {}, Return[{}, Module]];
                        tripleSet = AssociationMap[True &, realZeroTriples];
                        pairs = Select[realZeroTriples, KeyExistsQ[tripleSet, Reverse[#]] &];
                        If[pairs === {}, Return[{}, Module]];
                        edges = DeleteDuplicates[
                            UndirectedEdge[u[#[[1]], v], u[#[[3]], v]] & /@ pairs
                        ];
                        components = ConnectedComponents[Graph[edges]];
                        Flatten @ Map[
                            Function[comp,
                                With[{canonical = First @ Sort[comp]},
                                    Rule[#, canonical] & /@ DeleteCases[comp, canonical]
                                ]
                            ],
                            components
                        ]
                    ]
                ],
                verticesList
            ];

            altOptCond = altSwitch[j, u, scFn] @@@ activeTriples;
        ];

        mfgComplementarityData @ <|
            "SwitchingCosts"         -> switchingCosts,
            "ConsistentCosts"        -> consistCosts,
            "AltFlows"               -> altFlows,
            "AltTransitionFlows"     -> altTransitionFlows,
            "IneqSwitchingByVertex"  -> ineqSwitchingByVertex,
            "AltOptCond"             -> altOptCond,
            "ActiveTriples"          -> activeTriples,
            "ZeroSwitchUEqualities"  -> zeroSwitchUEqualities
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
        (* The Hamiltonian equation u[a,b]-u[b,a]+j[a,b]-j[b,a] = nrhs is enforced
           unconditionally. Dropping the former (j[a,b]-j[b,a]==0) || disjunct fixes
           a gap: on zero-flow edges that disjunct fired True, leaving u[a,b] and u[b,a]
           unconstrained. That prevented exit values at bypassed exit nodes from being
           pinned, producing parametric solutions with free value variables at unused
           exits. Without the disjunct, zero net flow forces u[a,b]=u[b,a]; the
           switching inequalities at adjacent nodes then propagate this to pin values
           at unused exits to their terminal cost. Safe for alpha≠1: nrhs=0 when j=0
           for any alpha, so the behaviour on zero-flow edges is unchanged. *)
        eqGeneral = MapThread[
            With[{eq = Equal[#1, If[alphaAtEdge[halfPairs[[#3]]] === 1, 0, #2]]},
                eq
            ] &,
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
         ruleEntryIn, ruleEntryOut, bm, km, vars},
        eqEntryIn            = systemData[sys, "EqEntryIn"];
        balanceGatheringFlows = systemData[sys, "BalanceGatheringFlows"];
        balanceSplittingFlows = systemData[sys, "BalanceSplittingFlows"];
        ruleEntryIn          = systemData[sys, "RuleEntryIn"];
        ruleEntryOut         = systemData[sys, "RuleEntryOut"];

        If[MissingQ[eqEntryIn],            eqEntryIn = True];
        If[MissingQ[balanceGatheringFlows], balanceGatheringFlows = {}];
        If[MissingQ[balanceSplittingFlows], balanceSplittingFlows = {}];
        If[MissingQ[ruleEntryIn],          ruleEntryIn = <||>];
        If[MissingQ[ruleEntryOut],         ruleEntryOut = <||>];

        kirchhoff = Join[
            If[Head[eqEntryIn] === List, eqEntryIn, {eqEntryIn}],
            (# == 0& /@ Join[balanceGatheringFlows, balanceSplittingFlows])
        ];
        kirchhoff = kirchhoff /. Join[Normal[ruleEntryIn], Normal[ruleEntryOut]];
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

(* When Alpha == 1, the generated system consists of linear equations,
   alternatives between linear equations, and linear inequalities; nothing is
   non-linear. *)
makeSystem[s_?scenarioQ] := makeSystem[s, makeSymbolicUnknowns[s]];

makeSystem[s_?scenarioQ, unk_?symbolicUnknownsQ] :=
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
        unknownAuxTriples = symbolicUnknownsData[unk, "AuxTriples"];
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
            (* Only the model fields accessed downstream via systemData are kept here.
               Vertices/Adjacency/Switching remain in scenarioData[s, "Model"]. *)
            KeyTake[model, {"Entries", "Exits"}],
            (* Only topology fields accessed downstream; AuxiliaryGraph/AuxEntryEdges/
               AuxExitEdges are not accessed via systemData. *)
            KeyTake[topology, {"Graph", "AuxEntryVertices", "AuxExitVertices",
                               "HalfPairs", "InAuxEntryPairs", "OutAuxExitPairs",
                               "Pairs", "AuxPairs"}],
            <|
                "Js"          -> symbolicUnknownsData[unk, "Js"],
                "Us"          -> symbolicUnknownsData[unk, "Us"],
                "Jts"         -> symbolicUnknownsData[unk, "Jts"],
                (* Keep modeling metadata available to solver-layer eligibility checks. *)
                "Hamiltonian" -> scenarioData[s, "Hamiltonian"],
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
