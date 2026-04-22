(*
   Scenario: typed scenario kernel for MFGraphs.

   Provides a typed wrapper around network-plus-parameters data so that
   scenarios can be constructed, validated, completed, and stored uniformly.

   Lifecycle:
     makeScenario[rawAssoc]  →  validates + completes + wraps in scenario[...] head
     validateScenario[s]     →  checks required keys; returns s or Failure[...]
     completeScenario[s]     →  fills Identity (hash), Benchmark defaults; returns s
     scenarioQ[x]            →  True iff x is a well-formed typed scenario

   Canonical top-level blocks inside scenario[<|...|>]:
     "Identity"    — name, version, contentHash
     "Model"       — raw network topology accepted by DataToEquations
     "Data"        — parameter substitution rules (e.g. {I1 -> 100, U1 -> 0})
     "Validation"  — consistency check results (populated after solving)
     "Benchmark"   — tier, timeout
     "Visualization" — optional plotting hints
     "Inheritance" — optional lineage/parent reference
*)

(* --- Public API declarations --- *)

scenario::usage =
"scenario[assoc] is the typed head for a MFGraphs scenario object. Use makeScenario \
to construct one; use ScenarioData to access keys.";

scenarioQ::usage =
"scenarioQ[x] returns True if x is a typed scenario[assoc_Association] object, \
False otherwise.";

makeScenario::usage =
"makeScenario[assoc] constructs a typed scenario from a raw association. The input \
must contain a \"Model\" key whose value is a network topology association accepted by \
DataToEquations (required keys: \"Vertices List\", \"Adjacency Matrix\", \
\"Entrance Vertices and Flows\", \"Exit Vertices and Terminal Costs\", \
\"Switching Costs\"). If \"Model\" contains a Wolfram Graph under key \"Graph\", \
missing \"Vertices List\" and/or \"Adjacency Matrix\" are derived automatically. \
Optional keys: \"Data\" (parameter substitution rules), \
\"Identity\" (name, version), \"Benchmark\" (tier, timeout), \"Visualization\", \
\"Inheritance\". Returns a scenario[...] object on success or Failure[...] on error.";

validateScenario::usage =
"validateScenario[s] checks that the scenario s has all required Model keys and that \
the Model value is an Association. Returns s unchanged on success, or \
Failure[\"ScenarioValidation\", <|\"Message\" -> msg, \"MissingKeys\" -> {...}|>] \
on failure.";

completeScenario::usage =
"completeScenario[s] fills in derived fields: sets \"contentHash\" in the Identity \
block (SHA256 of the canonical Model string), and supplies default Benchmark values \
(\"Tier\" -> \"core\", \"Timeout\" -> 300) when missing. Returns a new scenario object.";

ScenarioData::usage =
"ScenarioData[s, key] returns the value associated with key in the scenario s, or \
Missing[\"KeyAbsent\", key] if absent. ScenarioData[s] returns the underlying \
Association.";

(* Shared Topology Logic *)

BuildAuxiliaryTopology::usage =
"BuildAuxiliaryTopology[model] returns an association with the auxiliary graph, \
triples, pairs, and auxiliary vertex lists derived from the raw model.";

Begin["`Private`"];

(* Required keys that must appear inside the "Model" block. *)
$RequiredModelKeys = {
    "Vertices List",
    "Adjacency Matrix",
    "Entrance Vertices and Flows",
    "Exit Vertices and Terminal Costs",
    "Switching Costs"
};

(* Default benchmark settings applied by completeScenario. *)
$DefaultBenchmarkTier    = "core";
$DefaultBenchmarkTimeout = 300;

NormalizeScenarioModel[model_Association] :=
    Module[{graph, vertices, adjacency, normalized = model},
        graph = Lookup[model, "Graph", Missing["KeyAbsent", "Graph"]];
        If[!MatchQ[graph, _Graph], Return[model, Module]];

        vertices = Lookup[model, "Vertices List", Missing["KeyAbsent", "Vertices List"]];
        If[!ListQ[vertices],
            vertices = VertexList[graph]
        ];

        adjacency = Lookup[model, "Adjacency Matrix", Missing["KeyAbsent", "Adjacency Matrix"]];
        If[!ListQ[adjacency],
            adjacency = Normal @ AdjacencyMatrix[graph, vertices]
        ];

        If[!KeyExistsQ[normalized, "Vertices List"],
            normalized = Join[normalized, <|"Vertices List" -> vertices|>]
        ];
        If[!KeyExistsQ[normalized, "Adjacency Matrix"],
            normalized = Join[normalized, <|"Adjacency Matrix" -> adjacency|>]
        ];
        normalized
    ];

NormalizeScenarioModel[other_] := other;

(* --- Topology Helpers --- *)

BuildAuxiliaryTopology[model_Association] :=
    Module[{vertices, adjacency, entryFlows, exitCosts, graph, 
            entryVertices, exitVertices, auxEntryVertices, auxExitVertices,
            entryEdges, exitEdges, auxiliaryGraph, auxVerticesList, 
            auxTriples, edgeList, halfPairs, inAuxEntryPairs, 
            outAuxEntryPairs, inAuxExitPairs, outAuxExitPairs, pairs, auxPairs},
        
        vertices = Lookup[model, "Vertices List"];
        adjacency = Lookup[model, "Adjacency Matrix"];
        entryFlows = Lookup[model, "Entrance Vertices and Flows"];
        exitCosts = Lookup[model, "Exit Vertices and Terminal Costs"];
        
        If[!ListQ[vertices] || !ListQ[adjacency] || !ListQ[entryFlows] || !ListQ[exitCosts],
            Return[$Failed, Module]
        ];
        
        graph = Quiet @ Check[
            AdjacencyGraph[vertices, adjacency, DirectedEdges -> False],
            $Failed
        ];
        If[!MatchQ[graph, _Graph], Return[$Failed, Module]];

        entryVertices = First /@ entryFlows;
        exitVertices = First /@ exitCosts;
        
        auxEntryVertices = Symbol["en" <> ToString[#]] & /@ entryVertices;
        auxExitVertices = Symbol["ex" <> ToString[#]] & /@ exitVertices;
        
        entryEdges = MapThread[UndirectedEdge, {auxEntryVertices, entryVertices}];
        exitEdges = MapThread[UndirectedEdge, {exitVertices, auxExitVertices}];
        auxiliaryGraph = EdgeAdd[graph, Join[entryEdges, exitEdges]];
        
        auxVerticesList = VertexList[auxiliaryGraph];
        auxTriples = Flatten[
            Insert[#, 2] /@ Permutations[AdjacencyList[auxiliaryGraph, #], {2}] & /@ auxVerticesList,
            1
        ];

        edgeList = EdgeList[graph];
        halfPairs = List @@@ edgeList;
        inAuxEntryPairs = List @@@ entryEdges;
        outAuxExitPairs = List @@@ exitEdges;
        inAuxExitPairs = Reverse /@ outAuxExitPairs;
        outAuxEntryPairs = Reverse /@ inAuxEntryPairs;
        pairs = Join[halfPairs, Reverse /@ halfPairs];
        auxPairs = Join[inAuxEntryPairs, outAuxEntryPairs, inAuxExitPairs, outAuxExitPairs, pairs];

        <|
            "AuxiliaryGraph" -> auxiliaryGraph,
            "AuxTriples" -> auxTriples,
            "AuxPairs" -> auxPairs,
            "AuxEntryVertices" -> auxEntryVertices,
            "AuxExitVertices" -> auxExitVertices,
            "AuxEntryEdges" -> entryEdges,
            "AuxExitEdges" -> exitEdges
        |>
    ];

(* --- Type predicate --- *)

scenarioQ[scenario[_Association]] := True;
scenarioQ[_]                       := False;

(* --- Accessor --- *)

ScenarioData[scenario[assoc_Association]]           := assoc;
ScenarioData[scenario[assoc_Association], key_]     := Lookup[assoc, key, Missing["KeyAbsent", key]];

(* --- Validate --- *)

validateScenario[s_scenario] :=
    Module[{assoc, model, missing},
        assoc = ScenarioData[s];
        model = Lookup[assoc, "Model", Missing["KeyAbsent", "Model"]];
        Which[
            MissingQ[model],
                Failure["ScenarioValidation",
                    <|"Message" -> "\"Model\" key is absent.",
                      "MissingKeys" -> {"Model"}|>],
            !AssociationQ[model],
                Failure["ScenarioValidation",
                    <|"Message" -> "\"Model\" value must be an Association.",
                      "MissingKeys" -> {}|>],
            True,
                missing = Select[$RequiredModelKeys, !KeyExistsQ[model, #] &];
                If[missing === {},
                    s,
                    Failure["ScenarioValidation",
                        <|"Message" ->
                            "Missing required Model keys: " <> StringRiffle[missing, ", "],
                          "MissingKeys" -> missing|>]
                ]
        ]
    ];

(* Reject non-scenario inputs immediately. *)
validateScenario[x_] :=
    Failure["ScenarioValidation",
        <|"Message" -> "Input is not a scenario object.",
          "MissingKeys" -> {}|>];

(* --- Complete --- *)

CompleteSwitchingCosts[model_Association] :=
    Module[{topology, inputSC, scAssoc},
        topology = BuildAuxiliaryTopology[model];
        If[topology === $Failed, Return[model, Module]];
        
        inputSC = Lookup[model, "Switching Costs", {}];
        scAssoc = If[AssociationQ[inputSC],
            inputSC,
            AssociationThread[Most /@ inputSC, Last /@ inputSC]
        ];
        
        (* Explicitly complement with 0 for all possible transitions in the auxiliary graph *)
        Join[AssociationMap[0&, topology["AuxTriples"]], scAssoc]
    ];

completeScenario[s_scenario] :=
    Module[{assoc, identity, benchmark, model, hash, newAssoc, completedSC},
        assoc     = ScenarioData[s];
        model     = Lookup[assoc, "Model", <||>];
        identity  = Lookup[assoc, "Identity", <||>];
        benchmark = Lookup[assoc, "Benchmark", <||>];

        (* Ensure "Switching Costs" is explicitly completed before hashing. *)
        completedSC = CompleteSwitchingCosts[model];
        model = Join[model, <|"Switching Costs" -> completedSC|>];

        (* Compute content hash from the canonical string of the completed Model. *)
        hash = Hash[ToString[model, InputForm], "SHA256", "HexString"];

        (* Computed hash always wins: placed on the right so it overrides any user-supplied value. *)
        identity  = Join[identity, <|"contentHash" -> hash|>];
        benchmark = Join[
            <|"Tier" -> $DefaultBenchmarkTier, "Timeout" -> $DefaultBenchmarkTimeout|>,
            benchmark
        ];

        newAssoc = Join[assoc, <|
            "Model" -> model,
            "Identity" -> identity, 
            "Benchmark" -> benchmark
        |>];
        scenario[newAssoc]
    ];

completeScenario[x_] := x;   (* pass-through for non-scenario values *)

(* --- Constructor --- *)

makeScenario[rawAssoc_Association] :=
    Module[{normalizedAssoc, wrapped, validated, completed},
        normalizedAssoc = rawAssoc;
        If[AssociationQ[Lookup[rawAssoc, "Model", Missing["KeyAbsent", "Model"]]],
            normalizedAssoc = Join[
                rawAssoc,
                <|"Model" -> NormalizeScenarioModel[rawAssoc["Model"]]|>
            ]
        ];
        wrapped   = scenario[normalizedAssoc];
        validated = validateScenario[wrapped];
        If[FailureQ[validated],
            validated,
            completed = completeScenario[validated];
            completed
        ]
    ];

makeScenario[_] :=
    Failure["ScenarioValidation",
        <|"Message" -> "makeScenario requires an Association as input.",
          "MissingKeys" -> {}|>];

End[];
