(* Wolfram Language package *)
(*
   scenarioTools.wl — typed scenario kernel for MFGraphs.

   Provides a typed wrapper around network-plus-parameters data so that
   scenarios can be constructed, validated, completed, and stored uniformly.

   Lifecycle:
     makeScenario[rawAssoc]  →  validates + completes + wraps in scenario[...] head
     validateScenario[s]     →  checks required keys; returns s or Failure[...]
     completeScenario[s]     →  fills defaults and returns s
     scenarioQ[x]            →  True iff x is a well-formed typed scenario

   Canonical top-level blocks inside scenario[<|...|>]:
     "Identity"    — name, version
     "Model"       — raw network topology (keys: "Vertices", "Adjacency",
                     "Entries", "Exits", "Switching")
     "Validation"  — consistency check results
     "Benchmark"   — tier, timeout
     "Visualization" — optional plotting hints
     "Inheritance" — optional lineage/parent reference
     "Topology"    — derived auxiliary graph data (cached)
*)

BeginPackage["scenarioTools`", {"primitives`", "utilities`"}]

(* --- Public API declarations --- *)

scenario::usage =
"scenario[assoc] is the typed head for a MFGraphs scenario object. Use makeScenario \
to construct one; use scenarioData to access keys.";

scenarioQ::usage =
"scenarioQ[x] returns True if x is a typed scenario[assoc_Association] object, \
False otherwise.";

makeScenario::usage =
"makeScenario[assoc] constructs a typed scenario from a raw association. The input \
must contain a \"Model\" key whose value is a network topology association \
(required keys: \"Vertices\", \"Adjacency\", \"Entries\", \"Exits\", \"Switching\"). \
If \"Model\" contains a Wolfram Graph under key \"Graph\", missing \"Vertices\" \
and/or \"Adjacency\" are derived automatically. \
Optional keys: \
\"Hamiltonian\" (<|\"Alpha\" -> a, \"V\" -> v, \"G\" -> g, \
\"EdgeAlpha\" -> <|{u,v} -> a_uv, ...|>, \"EdgeV\" -> <|{u,v} -> v_uv, ...|>, \
\"EdgeG\" -> <|{u,v} -> g_uv, ...|>|>), \
\"Identity\" (name, version), \"Benchmark\" (tier, timeout), \"Visualization\", \
\"Inheritance\". Default Hamiltonian is Alpha=1 and V=-1 on all edges, with \
G[z]=-1/z (overridable globally and per edge). V/G/EdgeV/EdgeG are validated \
and preserved for future density and visualization work, but current system \
construction applies only Alpha/EdgeAlpha. Boundary values must be numeric; \
switching-cost values must be numeric or Infinity. Returns a scenario[...] object \
on success or Failure[...] on error.";

validateScenario::usage =
"validateScenario[s] checks that the scenario s has all required Model keys and that \
the Model value is an Association. Returns s unchanged on success, or \
Failure[\"ScenarioValidation\", <|\"Message\" -> msg, \"MissingKeys\" -> {...}|>] \
on failure.";

completeScenario::usage =
"completeScenario[s] fills in derived fields and supplies default Benchmark values \
(\"Tier\" -> \"core\", \"Timeout\" -> 300) when missing. If cached topology is \
absent, it warns and rebuilds topology from the Model before completing. Returns \
a new scenario object.";

scenarioData::usage =
"scenarioData[s, key] returns the value associated with key in the scenario s, or \
Missing[\"KeyAbsent\", key] if absent. scenarioData[s] returns the underlying \
Association.";

(* Shared Topology Logic *)

buildAuxiliaryTopology::usage =
"buildAuxiliaryTopology[model] returns an association with the auxiliary graph and \
metadata derived from the raw model.";

deriveAuxPairs::usage =
"deriveAuxPairs[topology] returns the list of all directed edge pairs {u, v} in the \
auxiliary graph (including entry/exit and reversed graph edges).";

buildAuxTriples::usage =
"buildAuxTriples[auxGraph] returns the list of all possible {v_in, v_mid, v_out} \
transitions (triples) in the graph.";

Begin["`Private`"]

scenarioFailure::usage =
"scenarioFailure[msg] builds a ScenarioValidation Failure. scenarioFailure[msg, keys] includes missing model keys.";

hamiltonianGTermQ::usage =
"hamiltonianGTermQ[value] returns True when value is a numeric Hamiltonian G term or a pure Function.";

buildAdjacencyFromGraph::usage =
"buildAdjacencyFromGraph[graph, vertices] returns graph's adjacency matrix ordered by vertices, or $Failed for unknown vertices.";

normalizeScenarioModel::usage =
"normalizeScenarioModel[model] fills Vertices and Adjacency from Graph when possible and normalizes the stored Graph.";

boundaryValuesNumericQ::usage =
"boundaryValuesNumericQ[model] returns True when Entries and Exits are vertex-value pairs with numeric values.";

switchingCostsNumericQ::usage =
"switchingCostsNumericQ[model] returns True when Switching has well-formed triples with numeric costs.";

normalizeSwitchingCosts::usage =
"normalizeSwitchingCosts[sc] converts list switching-cost rows to an Association and leaves Associations unchanged.";

integerVertexLabelsQ::usage =
"integerVertexLabelsQ[model] returns True when the model Vertices list contains only integer labels.";

modelDirectedEdgePairs::usage =
"modelDirectedEdgePairs[model] returns directed edge pairs derived from the model Vertices and Adjacency.";

normalizeHamiltonianSpec::usage =
"normalizeHamiltonianSpec[spec, model] expands Hamiltonian defaults and validates per-edge Hamiltonian overrides.";

buildAuxTriplesImpl::usage =
"buildAuxTriplesImpl[directedEdges] builds admissible auxiliary transition triples from directed edge pairs.";

completeSwitchingCosts::usage =
"completeSwitchingCosts[model, topology] fills absent switching costs with zero on the auxiliary topology.";

(* Required keys that must appear inside the "Model" block. *)
$RequiredModelKeys = {
    "Vertices",
    "Adjacency",
    "Entries",
    "Exits",
    "Switching"
};

(* Default benchmark settings applied by completeScenario. *)
$DefaultBenchmarkTier    = "core";
$DefaultBenchmarkTimeout = 300;
(* Default Hamiltonian:
   alpha = 1 on all edges, V = -1 on all edges, and g(z) = -1/z. *)
$DefaultHamiltonian = <|
    "Alpha" -> 1,
    "V" -> -1,
    "G" -> Function[z, -1 / z],
    "EdgeAlpha" -> <||>,
    "EdgeV" -> <||>,
    "EdgeG" -> <||>
|>;

scenarioFailure[msg_String] :=
    Failure["ScenarioValidation", <|"Message" -> msg, "MissingKeys" -> {}|>];
scenarioFailure[msg_String, missingKeys_List] :=
    Failure["ScenarioValidation", <|"Message" -> msg, "MissingKeys" -> missingKeys|>];

hamiltonianGTermQ[value_] :=
    NumericQ[value] || MatchQ[value, _Function];

buildAdjacencyFromGraph[graph_Graph, vertices_List] :=
    Module[{baseVertices, baseAdjacency, baseIndex, order},
        baseVertices = VertexList[graph];
        If[vertices === baseVertices,
            Return[AdjacencyMatrix[graph], Module]
        ];
        baseAdjacency = AdjacencyMatrix[graph];
        baseIndex = AssociationThread[baseVertices, Range[Length[baseVertices]]];
        order = Lookup[baseIndex, vertices, Missing["UnknownVertex"]];
        If[MemberQ[order, _Missing],
            $Failed,
            baseAdjacency[[order, order]]
        ]
    ];

normalizeScenarioModel[model_Association] :=
    Module[{graph, vertices, adjacency, normalized = model},
        graph = Lookup[model, "Graph", Missing["KeyAbsent", "Graph"]];
        If[!MatchQ[graph, _Graph], Return[model, Module]];

        vertices = Lookup[model, "Vertices", Missing["KeyAbsent", "Vertices"]];
        If[!ListQ[vertices],
            vertices = VertexList[graph]
        ];

        adjacency = Lookup[model, "Adjacency", Missing["KeyAbsent", "Adjacency"]];
        If[!(ListQ[adjacency] || Head[adjacency] === SparseArray),
            adjacency = buildAdjacencyFromGraph[graph, vertices]
        ];

        If[!KeyExistsQ[normalized, "Vertices"],
            normalized = Join[normalized, <|"Vertices" -> vertices|>]
        ];
        If[!KeyExistsQ[normalized, "Adjacency"],
            normalized = Join[normalized, <|"Adjacency" -> adjacency|>]
        ];
        normalized = Join[normalized, <|"Graph" -> AdjacencyGraph[vertices, Unitize[adjacency + Transpose[adjacency]], DirectedEdges -> False]|>];
        normalized
    ];

normalizeScenarioModel[other_] := other;

boundaryValuesNumericQ[model_Association] :=
    With[{entry = Lookup[model, "Entries", Missing[]],
          exit  = Lookup[model, "Exits", Missing[]]},
        ListQ[entry] && ListQ[exit] &&
        AllTrue[entry, MatchQ[#, {_, _}] &] &&
        AllTrue[exit,  MatchQ[#, {_, _}] &] &&
        AllTrue[Last /@ Join[entry, exit], NumericQ]
    ];

(* Boundary vertex membership is checked after integerVertexLabelsQ, so a separate
   IntegerQ test here would be redundant: every valid member of "Vertices" is
   already known to be an integer label. *)
boundaryVerticesPresentQ[model_Association] :=
    Module[{vertices, entryVertices, exitVertices},
        vertices = Lookup[model, "Vertices", Missing[]];
        If[!ListQ[vertices], Return[False, Module]];
        entryVertices = First /@ Lookup[model, "Entries", {}];
        exitVertices  = First /@ Lookup[model, "Exits", {}];
        AllTrue[Join[entryVertices, exitVertices], MemberQ[vertices, #] &]
    ];

(* Infinity is the explicit sentinel for blocked transitions. In WL this is the
   canonical positive infinity value; other directed infinities are intentionally
   not accepted as switching costs. *)
switchingCostValueQ[value_] :=
    NumericQ[value] || value === Infinity;

switchingCostsNumericQ[model_Association] :=
    Module[{switching},
        switching = Lookup[model, "Switching", Missing[]];
        Which[
            AssociationQ[switching],
                AllTrue[Keys[switching], ListQ[#] && Length[#] === 3 && AllTrue[#, IntegerQ] &] &&
                AllTrue[Values[switching], switchingCostValueQ],
            ListQ[switching],
                switching === {} ||
                (AllTrue[switching, ListQ[#] && Length[#] === 4 && AllTrue[Take[#, 3], IntegerQ] &] &&
                 AllTrue[Last /@ switching, switchingCostValueQ]),
            True,
                False
        ]
    ];

disjointBoundaryVerticesQ[model_Association] :=
    Module[{entryVertices, exitVertices},
        entryVertices = First /@ Lookup[model, "Entries", {}];
        exitVertices  = First /@ Lookup[model, "Exits", {}];
        Intersection[entryVertices, exitVertices] === {}
    ];

normalizeSwitchingCosts[sc_Association] := sc;
normalizeSwitchingCosts[sc_List] :=
    If[sc === {}, <||>, AssociationThread[Most /@ sc, Last /@ sc]];

integerVertexLabelsQ[model_Association] :=
    With[{vertices = Lookup[model, "Vertices", Missing[]]},
        ListQ[vertices] && AllTrue[vertices, IntegerQ]
    ];

modelDirectedEdgePairs[model_Association] :=
    Module[{vertices, adjacency},
        vertices  = Lookup[model, "Vertices", {}];
        adjacency = Lookup[model, "Adjacency", {}];
        If[vertices === {} || adjacency === {},
            {},
            List @@@ EdgeList[AdjacencyGraph[vertices, adjacency, DirectedEdges -> True]]
        ]
    ];

normalizeHamiltonianSpec[spec_, model_Association] :=
    Module[
        {ham, alphaDefault, vDefault, gDefault, edgeAlpha, edgeV, edgeG,
         edgePairs, validEdges, badEdges, badAlpha, badV, badG},
        ham          = If[AssociationQ[spec], spec, <||>];
        alphaDefault = Lookup[ham, "Alpha",     $DefaultHamiltonian["Alpha"]];
        vDefault     = Lookup[ham, "V",         $DefaultHamiltonian["V"]];
        gDefault     = Lookup[ham, "G",         $DefaultHamiltonian["G"]];
        edgeAlpha    = Lookup[ham, "EdgeAlpha", $DefaultHamiltonian["EdgeAlpha"]];
        edgeV        = Lookup[ham, "EdgeV",     $DefaultHamiltonian["EdgeV"]];
        edgeG        = Lookup[ham, "EdgeG",     $DefaultHamiltonian["EdgeG"]];

        If[!NumericQ[alphaDefault],
            Return[scenarioFailure["\"Hamiltonian\" default \"Alpha\" must be numeric."], Module]];
        If[!NumericQ[vDefault],
            Return[scenarioFailure["\"Hamiltonian\" default \"V\" must be numeric."], Module]];
        If[!hamiltonianGTermQ[gDefault],
            Return[scenarioFailure["\"Hamiltonian\" default \"G\" must be numeric or a pure function."], Module]];
        If[!AssociationQ[edgeAlpha],
            Return[scenarioFailure["\"Hamiltonian\" \"EdgeAlpha\" must be an Association."], Module]];
        If[!AssociationQ[edgeV],
            Return[scenarioFailure["\"Hamiltonian\" \"EdgeV\" must be an Association."], Module]];
        If[!AssociationQ[edgeG],
            Return[scenarioFailure["\"Hamiltonian\" \"EdgeG\" must be an Association."], Module]];

        edgePairs  = DeleteDuplicates[modelDirectedEdgePairs[model]];
        validEdges = AssociationThread[edgePairs, ConstantArray[True, Length[edgePairs]]];
        badEdges   = Select[
            DeleteDuplicates @ Join[Keys[edgeAlpha], Keys[edgeV], Keys[edgeG]],
            !MatchQ[#, {_, _}] || (!KeyExistsQ[validEdges, #] && !KeyExistsQ[validEdges, Reverse[#]]) &
        ];
        If[badEdges =!= {},
            Return[Failure["ScenarioValidation",
                <|"Message" -> "\"Hamiltonian\" edge-parameter keys must be valid edge pairs {u,v}.",
                  "MissingKeys" -> {}, "InvalidEdges" -> badEdges|>], Module]];

        badAlpha = Select[Normal[edgeAlpha], !NumericQ[Last[#]] &];
        badV     = Select[Normal[edgeV],     !NumericQ[Last[#]] &];
        badG     = Select[Normal[edgeG],     !hamiltonianGTermQ[Last[#]] &];
        If[badAlpha =!= {} || badV =!= {} || badG =!= {},
            Return[Failure["ScenarioValidation",
                <|"Message" -> "\"Hamiltonian\" edge-parameter values must be numeric (or Function for EdgeG).",
                  "MissingKeys" -> {}, "InvalidRules" -> Join[badAlpha, badV, badG]|>], Module]];

        Join[ham, <|
            "Alpha"     -> alphaDefault,
            "V"         -> vDefault,
            "G"         -> gDefault,
            "EdgeAlpha" -> edgeAlpha,
            "EdgeV"     -> edgeV,
            "EdgeG"     -> edgeG
        |>]
    ];

(* --- Topology Helpers --- *)

buildAuxTriplesImpl[directedEdges_List] :=
    Module[{incomingByVertex, outgoingByVertex, middleVertices, triples,
            auxEntryVertexQ, auxExitVertexQ},
        auxEntryVertexQ[v_] := StringQ[v] && StringStartsQ[v, "auxEntry"];
        auxExitVertexQ[v_]  := StringQ[v] && StringStartsQ[v, "auxExit"];
        incomingByVertex = GroupBy[directedEdges, Last -> First];
        outgoingByVertex = GroupBy[directedEdges, First -> Last];
        middleVertices   = Intersection[Keys[incomingByVertex], Keys[outgoingByVertex]];
        triples = DeleteDuplicates @ Flatten[
            Table[
                If[vIn =!= vOut, {vIn, vMid, vOut}, Nothing],
                {vMid, middleVertices},
                {vIn,  Lookup[incomingByVertex, vMid, {}]},
                {vOut, Lookup[outgoingByVertex, vMid, {}]}
            ],
            2
        ];
        Select[triples, !(auxExitVertexQ[First[#]] || auxEntryVertexQ[Last[#]]) &]
    ];

buildAuxTriples[auxGraph_Graph] :=
    Module[{edges},
        edges = List @@@ EdgeList[auxGraph];
        buildAuxTriplesImpl[Join[edges, Reverse[edges, {2}]]]
    ];

deriveAuxPairs[topology_Association] :=
    Lookup[topology, "AuxPairs",
        Module[{graph, halfPairs, inAuxEntryPairs, outAuxExitPairs, pairs},
            graph           = topology["Graph"];
            halfPairs       = List @@@ EdgeList[graph];
            inAuxEntryPairs = List @@@ topology["AuxEntryEdges"];
            outAuxExitPairs = List @@@ topology["AuxExitEdges"];
            pairs           = Join[halfPairs, Reverse /@ halfPairs];
            Join[inAuxEntryPairs, outAuxExitPairs, pairs]
        ]
    ];

buildAuxiliaryTopology[model_Association] :=
    Module[{vertices, adjacency, adjacencyForGraph, entryFlows, exitCosts, graph,
            entryVertices, exitVertices, auxEntryVertices, auxExitVertices,
            entryEdges, exitEdges, auxiliaryGraph, auxTriples,
            halfPairs, inAuxEntryPairs, outAuxExitPairs, pairs, auxPairs},

        vertices    = Lookup[model, "Vertices"];
        adjacency   = Lookup[model, "Adjacency"];
        entryFlows  = Lookup[model, "Entries"];
        exitCosts   = Lookup[model, "Exits"];

        If[
            !ListQ[vertices] ||
            !(Head[adjacency] === SparseArray || MatrixQ[adjacency]) ||
            !ListQ[entryFlows] ||
            !ListQ[exitCosts],
            Return[$Failed, Module]
        ];

        adjacencyForGraph = Unitize[adjacency + Transpose[adjacency]];

        graph = Check[
            AdjacencyGraph[vertices, adjacencyForGraph, DirectedEdges -> False],
            $Failed
        ];
        If[!MatchQ[graph, _Graph], Return[$Failed, Module]];

        entryVertices    = First /@ entryFlows;
        exitVertices     = First /@ exitCosts;
        auxEntryVertices = ("auxEntry" <> ToString[#]) & /@ entryVertices;
        auxExitVertices  = ("auxExit"  <> ToString[#]) & /@ exitVertices;

        entryEdges = MapThread[DirectedEdge, {auxEntryVertices, entryVertices}];
        exitEdges  = MapThread[DirectedEdge, {exitVertices, auxExitVertices}];

        halfPairs       = List @@@ EdgeList[graph];
        inAuxEntryPairs = List @@@ entryEdges;
        outAuxExitPairs = List @@@ exitEdges;
        pairs           = Join[halfPairs, Reverse[halfPairs, {2}]];
        auxPairs        = Join[inAuxEntryPairs, outAuxExitPairs, pairs];

        auxiliaryGraph = EdgeAdd[graph, Join[entryEdges, exitEdges]];
        auxTriples     = buildAuxTriplesImpl[auxPairs];

        <|
            "Graph"             -> graph,
            "AuxiliaryGraph"    -> auxiliaryGraph,
            "AuxEntryVertices"  -> auxEntryVertices,
            "AuxExitVertices"   -> auxExitVertices,
            "AuxEntryEdges"     -> entryEdges,
            "AuxExitEdges"      -> exitEdges,
            "AuxTriples"        -> auxTriples,
            "HalfPairs"         -> halfPairs,
            "InAuxEntryPairs"   -> inAuxEntryPairs,
            "OutAuxExitPairs"   -> outAuxExitPairs,
            "Pairs"             -> pairs,
            "AuxPairs"          -> auxPairs
        |>
    ];

(* --- Type predicate --- *)

scenarioQ[x_] := mfgTypedQ[x, scenario];

(* --- Accessor --- *)

scenarioData[s_] := mfgData[s];
scenarioData[s_, key_] := mfgData[s, key];

(* --- Validate --- *)

validateScenario[s_scenario] :=
    Module[{assoc, model, missing},
        assoc = scenarioData[s];
        model = Lookup[assoc, "Model", Missing["KeyAbsent", "Model"]];
        Which[
            MissingQ[model],
                scenarioFailure["\"Model\" key is absent.", {"Model"}],
            !AssociationQ[model],
                scenarioFailure["\"Model\" value must be an Association."],
            True,
                missing = Select[$RequiredModelKeys, !KeyExistsQ[model, #] &];
                If[missing === {},
                    s,
                    scenarioFailure[
                        "Missing required Model keys: " <> StringRiffle[missing, ", "],
                        missing]
                ]
        ]
    ];

validateScenario[x_] :=
    scenarioFailure["Input is not a scenario object."];

(* --- Complete --- *)

completeSwitchingCosts[model_Association, topology_Association] :=
    Module[{scAssoc, triples},
        scAssoc = Lookup[model, "Switching", <||>];
        triples = Lookup[topology, "AuxTriples", buildAuxTriples[topology["AuxiliaryGraph"]]];
        Join[AssociationMap[0&, triples], scAssoc]
    ];

completeScenario::notscenario =
    "completeScenario expected a scenario object; got `1`. Returning input unchanged.";
completeScenario::notopology =
    "completeScenario was called without cached topology; rebuilding topology from the Model before completing derived fields.";

completeScenario[s_scenario] :=
    Module[{assoc, identity, benchmark, model, hamiltonian, topology, completedSC},
        assoc       = scenarioData[s];
        model       = Lookup[assoc, "Model",       <||>];
        hamiltonian = Lookup[assoc, "Hamiltonian", $DefaultHamiltonian];
        identity    = Lookup[assoc, "Identity",    <||>];
        benchmark   = Lookup[assoc, "Benchmark",   <||>];
        topology    = Lookup[assoc, "Topology",    Missing["KeyAbsent", "Topology"]];

        If[MissingQ[topology] && AssociationQ[model],
            Message[completeScenario::notopology];
            model = Join[model, <|"Switching" -> normalizeSwitchingCosts[Lookup[model, "Switching", <||>]]|>];
            topology = buildAuxiliaryTopology[model]
        ];

        If[AssociationQ[topology],
            completedSC = completeSwitchingCosts[model, topology];
            model = Join[model, <|"Switching" -> completedSC|>]
        ];

        benchmark = Join[
            <|"Tier" -> $DefaultBenchmarkTier, "Timeout" -> $DefaultBenchmarkTimeout|>,
            benchmark
        ];

        scenario[Join[assoc, <|
            "Model"       -> model,
            "Hamiltonian" -> hamiltonian,
            "Identity"    -> identity,
            "Benchmark"   -> benchmark
        |>]]
    ];

completeScenario[x_] := (Message[completeScenario::notscenario, x]; x);

(* --- Constructor --- *)

makeScenario[rawAssoc_Association] :=
    Module[{rawModel, normalizedAssoc, validated, validatedAssoc, model, topology, hamiltonian},
        rawModel = Lookup[rawAssoc, "Model", Missing["KeyAbsent", "Model"]];
        normalizedAssoc = If[AssociationQ[rawModel],
            Join[rawAssoc, <|"Model" -> normalizeScenarioModel[rawModel]|>],
            rawAssoc
        ];

        validated = validateScenario[scenario[normalizedAssoc]];
        If[FailureQ[validated], Return[validated, Module]];

        validatedAssoc = scenarioData[validated];

        If[KeyExistsQ[validatedAssoc, "Data"],
            Return[scenarioFailure[
                "\"Data\" key is no longer supported. Provide numeric values directly in \"Model\"."],
                Module]
        ];

        model = validatedAssoc["Model"];
        If[!AssociationQ[model],
            Return[scenarioFailure["Model must remain an Association."], Module]];
        If[!integerVertexLabelsQ[model],
            Return[scenarioFailure["\"Vertices\" must contain only integers."], Module]];
        If[!boundaryValuesNumericQ[model],
            Return[scenarioFailure["Boundary values must be numeric."], Module]];
        If[!boundaryVerticesPresentQ[model],
            Return[scenarioFailure["Entry and exit vertices must belong to \"Vertices\"."], Module]];
        If[!switchingCostsNumericQ[model],
            Return[scenarioFailure["Switching cost values must be numeric."], Module]];
        If[!disjointBoundaryVerticesQ[model],
            Return[scenarioFailure["Entry and exit vertices must be disjoint."], Module]];

        model = Join[model, <|"Switching" -> normalizeSwitchingCosts[model["Switching"]]|>];

        hamiltonian = normalizeHamiltonianSpec[Lookup[validatedAssoc, "Hamiltonian", <||>], model];
        If[FailureQ[hamiltonian], Return[hamiltonian, Module]];

        topology = buildAuxiliaryTopology[model];
        If[!AssociationQ[topology],
            Return[scenarioFailure["Model topology is invalid or could not be constructed."], Module]];

        completeScenario[scenario[Join[
            validatedAssoc,
            <|"Model" -> model, "Topology" -> topology, "Hamiltonian" -> hamiltonian|>
        ]]]
    ];

makeScenario[_] :=
    scenarioFailure["makeScenario requires an Association as input."];

End[]

EndPackage[]
