(*
   Scenario: typed scenario kernel for MFGraphs.

   Provides a typed wrapper around network-plus-parameters data so that
   scenarios can be constructed, validated, completed, and stored uniformly.

   Lifecycle:
     makeScenario[rawAssoc]  →  validates + completes + wraps in scenario[...] head
     validateScenario[s]     →  checks required keys; returns s or Failure[...]
     completeScenario[s]     →  fills defaults and returns s
     scenarioQ[x]            →  True iff x is a well-formed typed scenario

   Canonical top-level blocks inside scenario[<|...|>]:
     "Identity"    — name, version
     "Model"       — raw network topology accepted by the core scenario/system kernels
     "Validation"  — consistency check results
     "Benchmark"   — tier, timeout
     "Visualization" — optional plotting hints
     "Inheritance" — optional lineage/parent reference
*)

BeginPackage["MFGraphs`"];

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
the core scenario/system kernels (required keys: \"Vertices List\", \"Adjacency Matrix\", \
\"Entrance Vertices and Flows\", \"Exit Vertices and Terminal Costs\", \
\"Switching Costs\"). If \"Model\" contains a Wolfram Graph under key \"Graph\", \
missing \"Vertices List\" and/or \"Adjacency Matrix\" are derived automatically. \
Optional keys: \
\"Hamiltonian\" (<|\"Alpha\" -> a, \"V\" -> v, \"G\" -> g, \
\"EdgeAlpha\" -> <|{u,v} -> a_uv, ...|>, \"EdgeV\" -> <|{u,v} -> v_uv, ...|>, \
\"EdgeG\" -> <|{u,v} -> g_uv, ...|>|>), \
\"Identity\" (name, version), \"Benchmark\" (tier, timeout), \"Visualization\", \
\"Inheritance\". Default Hamiltonian is Alpha=1 and V=0 on all edges, with \
G[z]=-1/z (overridable globally and per edge). Boundary and switching-cost values must be numeric. \
Returns a scenario[...] object on success or Failure[...] on error.";

validateScenario::usage =
"validateScenario[s] checks that the scenario s has all required Model keys and that \
the Model value is an Association. Returns s unchanged on success, or \
Failure[\"ScenarioValidation\", <|\"Message\" -> msg, \"MissingKeys\" -> {...}|>] \
on failure.";

completeScenario::usage =
"completeScenario[s] fills in derived fields and supplies default Benchmark values \
(\"Tier\" -> \"core\", \"Timeout\" -> 300) when missing. Returns a new scenario object.";

ScenarioData::usage =
"ScenarioData[s, key] returns the value associated with key in the scenario s, or \
Missing[\"KeyAbsent\", key] if absent. ScenarioData[s] returns the underlying \
Association.";

(* Shared Topology Logic *)

BuildAuxiliaryTopology::usage =
"BuildAuxiliaryTopology[model] returns an association with the auxiliary graph and \
metadata derived from the raw model.";

DeriveAuxPairs::usage =
"DeriveAuxPairs[topology] returns the list of all directed edge pairs {u, v} in the \
auxiliary graph (including entry/exit and reversed graph edges).";

BuildAuxTriples::usage =
"BuildAuxTriples[auxGraph] returns the list of all possible {v_in, v_mid, v_out} \
transitions (triples) in the graph.";

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
(* Default Hamiltonian:
   alpha = 1 on all edges, V = 0 on all edges, and g(z) = -1/z. *)
$DefaultHamiltonian = <|
    "Alpha" -> 1,
    "V" -> 0,
    "G" -> Function[z, -1 / z],
    "EdgeAlpha" -> <||>,
    "EdgeV" -> <||>,
    "EdgeG" -> <||>
|>;

HamiltonianGTermQ[value_] :=
    NumericQ[value] || MatchQ[value, _Function];

BuildAdjacencyFromGraph[graph_Graph, vertices_List] :=
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

NormalizeScenarioModel[model_Association] :=
    Module[{graph, vertices, adjacency, normalized = model},
        graph = Lookup[model, "Graph", Missing["KeyAbsent", "Graph"]];
        If[!MatchQ[graph, _Graph], Return[model, Module]];

        vertices = Lookup[model, "Vertices List", Missing["KeyAbsent", "Vertices List"]];
        If[!ListQ[vertices],
            vertices = VertexList[graph]
        ];

        adjacency = Lookup[model, "Adjacency Matrix", Missing["KeyAbsent", "Adjacency Matrix"]];
        If[!(ListQ[adjacency] || Head[adjacency] === SparseArray),
            adjacency = BuildAdjacencyFromGraph[graph, vertices]
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

BoundaryValuesNumericQ[model_Association] :=
    Module[{entry, exit, entryVals, exitVals},
        entry = Lookup[model, "Entrance Vertices and Flows", Missing["KeyAbsent", "Entrance Vertices and Flows"]];
        exit = Lookup[model, "Exit Vertices and Terminal Costs", Missing["KeyAbsent", "Exit Vertices and Terminal Costs"]];
        If[!ListQ[entry] || !ListQ[exit],
            Return[False, Module]
        ];
        If[!AllTrue[entry, MatchQ[#, {_, _}] &] || !AllTrue[exit, MatchQ[#, {_, _}] &],
            Return[False, Module]
        ];
        entryVals = Last /@ entry;
        exitVals = Last /@ exit;
        AllTrue[Join[entryVals, exitVals], NumericQ]
    ];

SwitchingCostsNumericQ[model_Association] :=
    Module[{switching},
        switching = Lookup[model, "Switching Costs", Missing["KeyAbsent", "Switching Costs"]];
        Which[
            AssociationQ[switching],
                AllTrue[Values[switching], NumericQ],
            ListQ[switching],
                AllTrue[switching, ListQ] && AllTrue[Last /@ switching, NumericQ],
            True,
                False
        ]
    ];

IntegerVertexLabelsQ[model_Association] :=
    Module[{vertices},
        vertices = Lookup[model, "Vertices List", Missing["KeyAbsent", "Vertices List"]];
        ListQ[vertices] && AllTrue[vertices, IntegerQ]
    ];

ModelDirectedEdgePairs[model_Association] :=
    Module[{vertices, adjacency},
        vertices = Lookup[model, "Vertices List", {}];
        adjacency = Lookup[model, "Adjacency Matrix", {}];
        If[vertices === {} || adjacency === {},
            {},
            List @@@ EdgeList[AdjacencyGraph[vertices, adjacency, DirectedEdges -> True]]
        ]
    ];

NormalizeHamiltonianSpec[spec_, model_Association] :=
    Module[
        {
            ham, alphaDefault, vDefault, gDefault, edgeAlpha, edgeV, edgeG,
            edgePairs, validEdges, badEdges, badValues
        },
        ham = If[AssociationQ[spec], spec, <||>];
        alphaDefault = Lookup[ham, "Alpha", $DefaultHamiltonian["Alpha"]];
        vDefault = Lookup[ham, "V", $DefaultHamiltonian["V"]];
        gDefault = Lookup[ham, "G", $DefaultHamiltonian["G"]];
        edgeAlpha = Lookup[ham, "EdgeAlpha", $DefaultHamiltonian["EdgeAlpha"]];
        edgeV = Lookup[ham, "EdgeV", $DefaultHamiltonian["EdgeV"]];
        edgeG = Lookup[ham, "EdgeG", $DefaultHamiltonian["EdgeG"]];
        If[!NumericQ[alphaDefault],
            Return[
                Failure["ScenarioValidation",
                    <|"Message" -> "\"Hamiltonian\" default \"Alpha\" must be numeric.",
                      "MissingKeys" -> {}|>
                ],
                Module
            ]
        ];
        If[!NumericQ[vDefault],
            Return[
                Failure["ScenarioValidation",
                    <|"Message" -> "\"Hamiltonian\" default \"V\" must be numeric.",
                      "MissingKeys" -> {}|>
                ],
                Module
            ]
        ];
        If[!HamiltonianGTermQ[gDefault],
            Return[
                Failure["ScenarioValidation",
                    <|"Message" -> "\"Hamiltonian\" default \"G\" must be numeric or a pure function.",
                      "MissingKeys" -> {}|>
                ],
                Module
            ]
        ];
        If[!AssociationQ[edgeAlpha],
            Return[
                Failure["ScenarioValidation",
                    <|"Message" -> "\"Hamiltonian\" \"EdgeAlpha\" must be an Association.",
                      "MissingKeys" -> {}|>
                ],
                Module
            ]
        ];
        If[!AssociationQ[edgeV],
            Return[
                Failure["ScenarioValidation",
                    <|"Message" -> "\"Hamiltonian\" \"EdgeV\" must be an Association.",
                      "MissingKeys" -> {}|>
                ],
                Module
            ]
        ];
        If[!AssociationQ[edgeG],
            Return[
                Failure["ScenarioValidation",
                    <|"Message" -> "\"Hamiltonian\" \"EdgeG\" must be an Association.",
                      "MissingKeys" -> {}|>
                ],
                Module
            ]
        ];
        edgePairs = DeleteDuplicates[ModelDirectedEdgePairs[model]];
        validEdges = AssociationThread[edgePairs, ConstantArray[True, Length[edgePairs]]];
        badEdges = Select[
            DeleteDuplicates @ Join[Keys[edgeAlpha], Keys[edgeV], Keys[edgeG]],
            !MatchQ[#, {_, _}] || (!KeyExistsQ[validEdges, #] && !KeyExistsQ[validEdges, Reverse[#]]) &
        ];
        If[badEdges =!= {},
            Return[
                Failure["ScenarioValidation",
                    <|"Message" -> "\"Hamiltonian\" edge-parameter keys must be valid edge pairs {u,v}.",
                      "MissingKeys" -> {},
                      "InvalidEdges" -> badEdges|>
                ],
                Module
            ]
        ];
        badValues = Select[
            Join[Normal[edgeAlpha], Normal[edgeV], Normal[edgeG]],
            !(
                NumericQ[Last[#]] ||
                (First[#] =!= {} && KeyExistsQ[edgeG, First[#]] && HamiltonianGTermQ[Last[#]])
            ) &
        ];
        If[badValues =!= {},
            Return[
                Failure["ScenarioValidation",
                    <|"Message" -> "\"Hamiltonian\" edge-parameter values must be numeric.",
                      "MissingKeys" -> {},
                      "InvalidRules" -> badValues|>
                ],
                Module
            ]
        ];
        Join[ham, <|
            "Alpha" -> alphaDefault,
            "V" -> vDefault,
            "G" -> gDefault,
            "EdgeAlpha" -> edgeAlpha,
            "EdgeV" -> edgeV,
            "EdgeG" -> edgeG
        |>]
    ];

(* --- Topology Helpers --- *)

iBuildAuxTriples[directedEdges_List] :=
    Module[{incomingByVertex, outgoingByVertex, middleVertices},
        incomingByVertex = GroupBy[directedEdges, Last -> First];
        outgoingByVertex = GroupBy[directedEdges, First -> Last];
        middleVertices = Intersection[Keys[incomingByVertex], Keys[outgoingByVertex]];
        DeleteDuplicates @ Flatten[
            Table[
                If[vIn =!= vOut, {vIn, vMid, vOut}, Nothing],
                {vMid, middleVertices},
                {vIn, Lookup[incomingByVertex, vMid, {}]},
                {vOut, Lookup[outgoingByVertex, vMid, {}]}
            ],
            2
        ]
    ];

BuildAuxTriples[auxGraph_Graph] :=
    Module[{edges},
        edges = List @@@ EdgeList[auxGraph];
        iBuildAuxTriples[Join[edges, Reverse[edges, {2}]]]
    ];

DeriveAuxPairs[topology_Association] :=
    Lookup[topology, "AuxPairs",
        Module[{graph, halfPairs, inAuxEntryPairs, outAuxExitPairs, inAuxExitPairs,
                outAuxEntryPairs, pairs},
            graph = topology["Graph"];
            halfPairs = List @@@ EdgeList[graph];
            inAuxEntryPairs = List @@@ topology["AuxEntryEdges"];
            outAuxExitPairs = List @@@ topology["AuxExitEdges"];
            inAuxExitPairs = Reverse /@ outAuxExitPairs;
            outAuxEntryPairs = Reverse /@ inAuxEntryPairs;
            pairs = Join[halfPairs, Reverse /@ halfPairs];
            Join[inAuxEntryPairs, outAuxEntryPairs, inAuxExitPairs, outAuxExitPairs, pairs]
        ]
    ];

BuildAuxiliaryTopology[model_Association] :=
    Module[{vertices, adjacency, adjacencyForGraph, entryFlows, exitCosts, graph, 
            entryVertices, exitVertices, auxEntryVertices, auxExitVertices,
            entryEdges, exitEdges, auxiliaryGraph, auxTriples,
            halfPairs, inAuxEntryPairs, outAuxExitPairs, inAuxExitPairs,
            outAuxEntryPairs, pairs, auxPairs},
        
        vertices = Lookup[model, "Vertices List"];
        adjacency = Lookup[model, "Adjacency Matrix"];
        entryFlows = Lookup[model, "Entrance Vertices and Flows"];
        exitCosts = Lookup[model, "Exit Vertices and Terminal Costs"];
        
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

        entryVertices = First /@ entryFlows;
        exitVertices = First /@ exitCosts;
        
        auxEntryVertices = Symbol["en" <> ToString[#]] & /@ entryVertices;
        auxExitVertices = Symbol["ex" <> ToString[#]] & /@ exitVertices;
        
        entryEdges = MapThread[DirectedEdge, {auxEntryVertices, entryVertices}];
        exitEdges = MapThread[DirectedEdge, {exitVertices, auxExitVertices}];
        
        halfPairs = List @@@ EdgeList[graph];
        inAuxEntryPairs = List @@@ entryEdges;
        outAuxExitPairs = List @@@ exitEdges;
        inAuxExitPairs = Reverse[outAuxExitPairs, {2}];
        outAuxEntryPairs = Reverse[inAuxEntryPairs, {2}];
        pairs = Join[halfPairs, Reverse[halfPairs, {2}]];
        
        auxPairs = Join[inAuxEntryPairs, outAuxEntryPairs, inAuxExitPairs, outAuxExitPairs, pairs];
        
        auxiliaryGraph = EdgeAdd[graph, Join[entryEdges, exitEdges]];
        auxTriples = iBuildAuxTriples[auxPairs];

        <|
            "Graph" -> graph,
            "AuxiliaryGraph" -> auxiliaryGraph,
            "AuxEntryVertices" -> auxEntryVertices,
            "AuxExitVertices" -> auxExitVertices,
            "AuxEntryEdges" -> entryEdges,
            "AuxExitEdges" -> exitEdges,
            "AuxTriples" -> auxTriples,
            "HalfPairs" -> halfPairs,
            "InAuxEntryPairs" -> inAuxEntryPairs,
            "OutAuxExitPairs" -> outAuxExitPairs,
            "InAuxExitPairs" -> inAuxExitPairs,
            "OutAuxEntryPairs" -> outAuxEntryPairs,
            "Pairs" -> pairs,
            "AuxPairs" -> auxPairs
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

CompleteSwitchingCosts[model_Association, topology_Association] :=
    Module[{inputSC, scAssoc, triples},
        inputSC = Lookup[model, "Switching Costs", {}];
        scAssoc = If[AssociationQ[inputSC],
            inputSC,
            AssociationThread[Most /@ inputSC, Last /@ inputSC]
        ];
        
        (* Derive triples for completion *)
        triples = Lookup[topology, "AuxTriples", BuildAuxTriples[topology["AuxiliaryGraph"]]];
        
        (* Explicitly complement with 0 for all possible transitions in the auxiliary graph *)
        Join[AssociationMap[0&, triples], scAssoc]
    ];

completeScenario[s_scenario] :=
    Module[{assoc, identity, benchmark, model, hamiltonian, newAssoc, completedSC, topology},
        assoc     = ScenarioData[s];
        model     = Lookup[assoc, "Model", <||>];
        hamiltonian = Lookup[assoc, "Hamiltonian", $DefaultHamiltonian];
        identity  = Lookup[assoc, "Identity", <||>];
        benchmark = Lookup[assoc, "Benchmark", <||>];
        topology  = Lookup[assoc, "Topology", Missing["KeyAbsent", "Topology"]];

        (* Ensure "Switching Costs" is explicitly completed before hashing. *)
        If[!MissingQ[topology],
            completedSC = CompleteSwitchingCosts[model, topology];
            model = Join[model, <|"Switching Costs" -> completedSC|>]
        ];

        benchmark = Join[
            <|"Tier" -> $DefaultBenchmarkTier, "Timeout" -> $DefaultBenchmarkTimeout|>,
            benchmark
        ];

        newAssoc = Join[assoc, <|
            "Model" -> model,
            "Hamiltonian" -> hamiltonian,
            "Identity" -> identity, 
            "Benchmark" -> benchmark
        |>];
        scenario[newAssoc]
    ];

completeScenario[x_] := x;   (* pass-through for non-scenario values *)

(* --- Constructor --- *)

makeScenario[rawAssoc_Association] :=
    Module[{normalizedAssoc, wrapped, validated, validatedAssoc, completed, model, topology, hamiltonian},
        normalizedAssoc = rawAssoc;
        model = Lookup[rawAssoc, "Model", Missing["KeyAbsent", "Model"]];
        
        If[AssociationQ[model],
            normalizedAssoc = Join[
                rawAssoc,
                <|"Model" -> NormalizeScenarioModel[model]|>
            ]
        ];
        
        wrapped   = scenario[normalizedAssoc];
        validated = validateScenario[wrapped];
        If[FailureQ[validated],
            validated,
            validatedAssoc = ScenarioData[validated];
            If[KeyExistsQ[validatedAssoc, "Data"],
                Return[
                    Failure["ScenarioValidation",
                        <|"Message" -> "\"Data\" key is no longer supported. Provide numeric values directly in \"Model\".",
                          "MissingKeys" -> {}|>
                    ],
                    Module
                ]
            ];
            model = validatedAssoc["Model"];
            If[!AssociationQ[model],
                Return[
                    Failure["ScenarioValidation",
                        <|"Message" -> "Model must remain an Association.",
                          "MissingKeys" -> {}|>
                    ],
                    Module
                ]
            ];
            If[!IntegerVertexLabelsQ[model],
                Return[
                    Failure["ScenarioValidation",
                        <|"Message" -> "\"Vertices List\" must contain only integers.",
                          "MissingKeys" -> {}|>
                    ],
                    Module
                ]
            ];
            If[!BoundaryValuesNumericQ[model],
                Return[
                    Failure["ScenarioValidation",
                        <|"Message" -> "Boundary values must be numeric.",
                          "MissingKeys" -> {}|>
                    ],
                    Module
                ]
            ];
            If[!SwitchingCostsNumericQ[model],
                Return[
                    Failure["ScenarioValidation",
                        <|"Message" -> "Switching cost values must be numeric.",
                          "MissingKeys" -> {}|>
                    ],
                    Module
                ]
            ];
            hamiltonian = NormalizeHamiltonianSpec[
                Lookup[validatedAssoc, "Hamiltonian", <||>],
                model
            ];
            If[FailureQ[hamiltonian],
                Return[hamiltonian, Module]
            ];
            topology = BuildAuxiliaryTopology[model];
            If[!AssociationQ[topology],
                Return[
                    Failure["ScenarioValidation",
                        <|"Message" -> "Model topology is invalid or could not be constructed.",
                          "MissingKeys" -> {}|>
                    ],
                    Module
                ]
            ];
            completed = completeScenario[scenario[Join[
                validatedAssoc,
                <|"Model" -> model, "Topology" -> topology, "Hamiltonian" -> hamiltonian|>
            ]]];
            completed
        ]
    ];

makeScenario[_] :=
    Failure["ScenarioValidation",
        <|"Message" -> "makeScenario requires an Association as input.",
          "MissingKeys" -> {}|>];

End[];

EndPackage[];
