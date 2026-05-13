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

boundaryPreservingAutomorphisms::usage = "boundaryPreservingAutomorphisms[scenario] returns the list of graph automorphisms of the scenario's underlying graph that also preserve the (vertex,value) multisets of Entries and Exits. Empty list if no nontrivial symmetry survives.";

addSymmetryEqualities::usage = "addSymmetryEqualities[sys, scenario] returns a new mfgSystem whose EqGeneral has the orbit equalities j[a,b] == j[sigma(a), sigma(b)] and u[a,b] == u[sigma(a), sigma(b)] appended for every boundary-preserving automorphism sigma. Symbolic, exact; no-op when the symmetry group is trivial.";

addOracleEqualities::usage = "addOracleEqualities[sys, oracleResult] returns a new mfgSystem with j[e] == 0 appended to EqGeneral for every e in oracleResult[\"Inactive\"]. Designed to consume the Association returned by numericOracleClassify (in the optional numericOracle subpackage); this function itself is symbolic-only. No-op when Inactive is empty or when oracleResult[\"Converged\"] is False. Runs a FindInstance safety check to verify the pruned system is still feasible; if not, returns sys unchanged.";

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
            If[zeroSwitchUEqualities =!= {},
                altOptCond = DeleteCases[altOptCond /. zeroSwitchUEqualities, True]
            ];
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
    Module[{halfPairs, signedFlows, nlhs, eqGeneral},
        (* The Hamiltonian equation u[a,b]-u[b,a]+j[a,b]-j[b,a] == 0 is enforced
           unconditionally. The (j[a,b]-j[b,a]==0) || ... disjunct must NOT be
           reintroduced: on zero-flow edges it fires True, leaving u[a,b]/u[b,a]
           unconstrained and producing parametric solutions with free value
           variables at unused exits. Without it, zero net flow forces
           u[a,b]==u[b,a], and switching inequalities propagate the value pin.
           Critical-congestion-only: under withCriticalCongestionGuard the per-edge
           alpha branching collapses to nlhs == 0; revisit if Alpha != 1 is ever
           supported. *)
        halfPairs   = topology["HalfPairs"];
        signedFlows = First[flowData]["SignedFlows"];
        nlhs = Flatten[u @@ # - u @@ Reverse @ # + signedFlows[#]& /@ halfPairs];
        eqGeneral = (# == 0) & /@ nlhs;

        mfgHamiltonianData @ <|
            "Nlhs"      -> nlhs,
            "EqGeneral" -> eqGeneral
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

        kirchhoff = Join[
            eqEntryIn,
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

(* --- Constructors --- *)

(* The generated system consists of linear equations, alternatives between
   linear equations, and linear inequalities; nothing is non-linear. Linearity
   is preserved by two conventions: (i) switching costs are scalars or affine
   Function[x, a x + b], so cost[j[r,i,w]] stays linear in j; (ii) non-linear
   Hamiltonian terms (e.g. m^alpha for Alpha != 1) are represented by cpf
   placeholder symbols, with the non-linear relationship deferred to a
   separate closure equation outside this kernel. *)
makeSystem[s_?scenarioQ] := makeSystem[s, makeSymbolicUnknowns[s]];

makeSystem[s_?scenarioQ, unk_?symbolicUnknownsQ] :=
    Module[{model, topology, graph, auxiliaryGraph, auxEdges, edges,
         auxVertices, auxTriples, unknownAuxTriples, boundaryData, flowData, compData, hamData,
         finalAssoc},

        model    = scenarioData[s, "Model"];
        topology = scenarioData[s, "Topology"];

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

(* ---- Boundary-aware symmetry folding ---- *)

(* Apply a Cycles permutation to a single vertex label. Integer vertices
   permute directly. Aux identifiers ("auxEntry<n>" / "auxExit<n>") are
   strings carrying the underlying vertex in their name — extract n,
   permute, reassemble so swapping vertex i with vertex j also swaps the
   matching aux. *)
applyPermVertex[v_Integer, perm_Cycles] := PermutationReplace[v, perm];
applyPermVertex[v_String, perm_Cycles] :=
    Module[{digits, prefix, n, permed},
        digits = StringCases[v, DigitCharacter ..];
        If[digits === {}, Return[v, Module]];
        prefix = StringDrop[v, -StringLength[Last[digits]]];
        n      = ToExpression[Last[digits]];
        permed = PermutationReplace[n, perm];
        If[permed === n, v, prefix <> ToString[permed]]
    ];
applyPermVertex[v_, _Cycles] := v;

(* Apply perm to a flow / value variable. j[a,b], j[a,b,c], u[a,b]
   permute their integer arguments componentwise. *)
applyPermVar[j[a_, b_],    perm_Cycles] := j[applyPermVertex[a, perm], applyPermVertex[b, perm]];
applyPermVar[j[a_, b_, c_], perm_Cycles] := j[applyPermVertex[a, perm], applyPermVertex[b, perm], applyPermVertex[c, perm]];
applyPermVar[u[a_, b_],    perm_Cycles] := u[applyPermVertex[a, perm], applyPermVertex[b, perm]];

(* True iff perm preserves the (vertex,value) multiset. Integer vertices
   get permuted; non-integers (aux symbols) stay put because applyPermVertex
   is the identity on them. *)
boundaryPreservedByQ[pairs_List, perm_Cycles] :=
    Sort[pairs] === Sort[{applyPermVertex[#[[1]], perm], #[[2]]} & /@ pairs];

boundaryPreservingAutomorphisms[s_?scenarioQ] :=
    Module[{graph, group, generators, entries, exits, model, valid},
        model = scenarioData[s, "Model"];
        graph = Lookup[scenarioData[s, "Topology"], "Graph", Null];
        If[graph === Null, Return[{}, Module]];
        entries = Lookup[model, "Entries", {}];
        exits   = Lookup[model, "Exits",   {}];
        group = Quiet @ TimeConstrained[GraphAutomorphismGroup[graph], 30, Null];
        If[group === Null || group === $Aborted, Return[{}, Module]];
        generators = If[Head[group] === PermutationGroup, GroupGenerators[group], {}];
        valid = Select[generators,
            boundaryPreservedByQ[entries, #] && boundaryPreservedByQ[exits, #] &];
        valid
    ];

(* Append j[e] == 0 to EqGeneral for every Inactive variable in the oracle
   result. Symbolic-only — input is a list of symbolic flow variables, output
   is a new mfgSystem with strictly more equalities. Conservative on three
   axes: empty Inactive list -> no-op; "Converged" -> False -> no-op; the
   pruned system fails its FindInstance safety check -> no-op. *)
addOracleEqualities[sys_?mfgSystemQ, oracle_Association] :=
    Module[{inactive, eqs, sysAssoc, hamRecord, hamAssoc, newEqGeneral, newHam,
            newSys, allVars, feasibilityProbe},
        If[! TrueQ[Lookup[oracle, "Converged", False]], Return[sys, Module]];
        inactive = Lookup[oracle, "Inactive", {}];
        If[inactive === {}, Return[sys, Module]];
        eqs = (# == 0) & /@ inactive;
        sysAssoc  = mfgData[sys];
        hamRecord = Lookup[sysAssoc, "HamiltonianData", Null];
        If[! mfgHamiltonianDataQ[hamRecord], Return[sys, Module]];
        hamAssoc  = mfgData[hamRecord];
        newEqGeneral = Join[Lookup[hamAssoc, "EqGeneral", {}], eqs];
        newHam = mfgHamiltonianData[<|hamAssoc, "EqGeneral" -> newEqGeneral|>];
        newSys = mfgSystem[<|sysAssoc, "HamiltonianData" -> newHam|>];
        (* Safety check: does the pruned linear part still admit a feasible
           point? If FindInstance returns {} the oracle over-pruned -- bail
           out and return the original system. Skip the check if the system
           lacks the expected linear-program buckets. *)
        allVars = Join[
            Lookup[sysAssoc, "Js",  {}],
            Lookup[sysAssoc, "Jts", {}],
            Lookup[sysAssoc, "Us",  {}]
        ];
        feasibilityProbe = TimeConstrained[
            Quiet @ FindInstance[
                And @@ Join[
                    Lookup[sysAssoc, "EqEntryIn", {}],
                    {Lookup[sysAssoc, "EqBalanceSplittingFlows", True],
                     Lookup[sysAssoc, "EqBalanceGatheringFlows", True]},
                    Lookup[sysAssoc, "EqGeneral", {}],
                    newEqGeneral,
                    Lookup[sysAssoc, "IneqJs",  {}],
                    Lookup[sysAssoc, "IneqJts", {}]
                ],
                allVars, Reals, 1
            ],
            5.0,
            $Aborted
        ];
        Which[
            feasibilityProbe === $Aborted,           newSys,  (* trust it *)
            MatchQ[feasibilityProbe, {{___Rule}, ___}], newSys,
            True,                                     sys     (* over-pruned *)
        ]
    ];

addSymmetryEqualities[sys_?mfgSystemQ, s_?scenarioQ] :=
    Module[{perms, js, jts, us, eqs, hamRecord, hamAssoc, newEqGeneral, newHam, sysAssoc},
        perms = boundaryPreservingAutomorphisms[s];
        If[perms === {}, Return[sys, Module]];
        js  = systemData[sys, "Js"];
        jts = systemData[sys, "Jts"];
        us  = systemData[sys, "Us"];
        eqs = Flatten @ Table[
            Join[
                # == applyPermVar[#, perm] & /@ js,
                # == applyPermVar[#, perm] & /@ jts,
                # == applyPermVar[#, perm] & /@ us
            ],
            {perm, perms}
        ];
        eqs = DeleteDuplicates @ DeleteCases[eqs, lhs_ == rhs_ /; lhs === rhs];
        sysAssoc  = mfgData[sys];
        hamRecord = Lookup[sysAssoc, "HamiltonianData", Null];
        If[! mfgHamiltonianDataQ[hamRecord], Return[sys, Module]];
        hamAssoc  = mfgData[hamRecord];
        newEqGeneral = Join[Lookup[hamAssoc, "EqGeneral", {}], eqs];
        newHam = mfgHamiltonianData[<|hamAssoc, "EqGeneral" -> newEqGeneral|>];
        mfgSystem[<|sysAssoc, "HamiltonianData" -> newHam|>]
    ];

End[]

EndPackage[]
