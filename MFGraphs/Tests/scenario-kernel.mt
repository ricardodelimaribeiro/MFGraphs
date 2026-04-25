(* Smoke tests for the scenario kernel (makeScenario, validateScenario, completeScenario, scenarioQ). *)

(* Inline test fixture — case 12 (4-vertex attraction network) with concrete values. *)
$testModel12 = <|
    "Vertices List"                    -> {1, 2, 3, 4},
    "Adjacency Matrix"                 -> {{0,1,1,0},{0,0,1,1},{0,0,0,1},{0,0,0,0}},
    "Entrance Vertices and Flows"      -> {{1, 100}},
    "Exit Vertices and Terminal Costs" -> {{4, 0}},
    "Switching Costs"                  -> {}
|>;

(* Test: public symbols exist in MFGraphs` context *)
Test[
    NameQ["MFGraphs`scenario"] &&
    NameQ["MFGraphs`scenarioQ"] &&
    NameQ["MFGraphs`makeScenario"] &&
    NameQ["MFGraphs`validateScenario"] &&
    NameQ["MFGraphs`completeScenario"] &&
    NameQ["MFGraphs`ScenarioData"]
    ,
    True
    ,
    TestID -> "Scenario kernel: public symbols exist in MFGraphs context"
]

(* Test: notebook ergonomics — unqualified symbols resolve after Needs *)
Test[
    Module[{names, globalNames},
        names = {
            "scenario",
            "scenarioQ",
            "makeScenario",
            "validateScenario",
            "completeScenario",
            "ScenarioData"
        };
        globalNames = ("Global`" <> #) & /@ names;
        Scan[
            Function[s,
                If[NameQ[s], Remove[s]]
            ],
            globalNames
        ];
        Needs["MFGraphs`"];
        And @@ (NameQ /@ names)
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: unqualified symbols are available after Needs"
]
(* Test: makeScenario returns a typed scenario for valid input *)
Test[
    Module[{data, s},
        data = $testModel12;
        s = makeScenario[<|"Model" -> data|>];
        scenarioQ[s]
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: makeScenario returns typed scenario"
]

(* Test: Benchmark block defaults are filled *)
Test[
    Module[{data, s, bench},
        data = $testModel12;
        s = makeScenario[<|"Model" -> data|>];
        bench = ScenarioData[s, "Benchmark"];
        bench["Tier"] === "core" && bench["Timeout"] === 300
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: default Benchmark tier and timeout"
]

(* Test: Hamiltonian block defaults are filled *)
Test[
    Module[{data, s, h},
        data = $testModel12;
        s = makeScenario[<|"Model" -> data|>];
        h = ScenarioData[s, "Hamiltonian"];
        AssociationQ[h] &&
        h["Alpha"] === 1 &&
        h["V"] === 0 &&
        MatchQ[h["G"], _Function] &&
        h["G"][2] === -1/2 &&
        h["EdgeAlpha"] === <||> &&
        h["EdgeV"] === <||> &&
        h["EdgeG"] === <||>
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: default Hamiltonian settings are present"
]

(* Test: user-supplied per-edge Hamiltonian parameters are preserved *)
Test[
    Module[{model, s, h},
        model = <|
            "Vertices List" -> {1, 2, 3},
            "Adjacency Matrix" -> {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}},
            "Entrance Vertices and Flows" -> {{1, 10}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs" -> <||>
        |>;
        s = makeScenario[<|
            "Model" -> model,
            "Hamiltonian" -> <|
                "Alpha" -> 1,
                "V" -> 3,
                "G" -> 4,
                "EdgeAlpha" -> <|{1, 2} -> 2|>,
                "EdgeV" -> <|{1, 2} -> 5|>,
                "EdgeG" -> <|{1, 2} -> 6|>
            |>
        |>];
        h = ScenarioData[s, "Hamiltonian"];
        scenarioQ[s] &&
        h["EdgeAlpha"][{1, 2}] === 2 &&
        h["EdgeV"][{1, 2}] === 5 &&
        h["EdgeG"][{1, 2}] === 6 &&
        h["V"] === 3 &&
        h["G"] === 4
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: per-edge Hamiltonian parameters are preserved"
]

(* Test: function-valued per-edge G is preserved *)
Test[
    Module[{model, s, h, gfun},
        model = <|
            "Vertices List" -> {1, 2, 3},
            "Adjacency Matrix" -> {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}},
            "Entrance Vertices and Flows" -> {{1, 10}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs" -> <||>
        |>;
        gfun = Function[z, -2/z];
        s = makeScenario[<|
            "Model" -> model,
            "Hamiltonian" -> <|"EdgeG" -> <|{1, 2} -> gfun|>|>
        |>];
        h = ScenarioData[s, "Hamiltonian"];
        scenarioQ[s] &&
        MatchQ[h["EdgeG"][{1, 2}], _Function] &&
        h["EdgeG"][{1, 2}][2] === -1
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: function-valued EdgeG is accepted"
]

(* Test: user-supplied Benchmark values override defaults *)
Test[
    Module[{data, s, bench},
        data = $testModel12;
        s = makeScenario[<|"Model" -> data, "Benchmark" -> <|"Tier" -> "stress", "Timeout" -> 900|>|>];
        bench = ScenarioData[s, "Benchmark"];
        bench["Tier"] === "stress" && bench["Timeout"] === 900
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: user Benchmark values are preserved"
]

(* Test: missing required Model key produces Failure *)
Test[
    Module[{incomplete, result},
        incomplete = <|
            "Vertices List" -> {1, 2},
            "Adjacency Matrix" -> {{0, 1}, {0, 0}},
            "Entrance Vertices and Flows" -> {{1, 100}},
            "Exit Vertices and Terminal Costs" -> {{2, 0}}
            (* "Switching Costs" intentionally absent *)
        |>;
        result = makeScenario[<|"Model" -> incomplete|>];
        FailureQ[result] && result["Tag"] === "ScenarioValidation"
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: missing Model key yields Failure"
]

(* Test: absent Model key produces Failure *)
Test[
    Module[{result},
        result = makeScenario[<|"Identity" -> <|"name" -> "no model"|>|>];
        FailureQ[result]
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: absent Model key yields Failure"
]

(* Test: non-scenario input to validateScenario yields Failure *)
Test[
    FailureQ[validateScenario["not a scenario"]]
    ,
    True
    ,
    TestID -> "Scenario kernel: non-scenario input rejected by validateScenario"
]

(* Test: ScenarioData accessor returns underlying association *)
Test[
    Module[{data, s, underlying},
        data = $testModel12;
        s = makeScenario[<|"Model" -> data|>];
        underlying = ScenarioData[s];
        AssociationQ[underlying] && KeyExistsQ[underlying, "Model"]
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: ScenarioData returns association"
]

(* Test: partial Identity override — user name is preserved *)
Test[
    Module[{data, s, identity},
        data = $testModel12;
        s = makeScenario[<|"Model" -> data, "Identity" -> <|"name" -> "my-scenario"|>|>];
        identity = ScenarioData[s, "Identity"];
        identity["name"] === "my-scenario"
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: partial Identity override preserves name"
]

(* Test: sparse adjacency matrix is accepted and topology is materialized *)
Test[
    Module[{model, s, topology, normalizedModel},
        model = <|
            "Vertices List" -> {1, 2},
            "Adjacency Matrix" -> SparseArray[{{1, 2} -> 1, {2, 1} -> 1}, {2, 2}],
            "Entrance Vertices and Flows" -> {{1, 10}},
            "Exit Vertices and Terminal Costs" -> {{2, 0}},
            "Switching Costs" -> <||>
        |>;
        s = makeScenario[<|"Model" -> model|>];
        topology = ScenarioData[s, "Topology"];
        normalizedModel = ScenarioData[s, "Model"];
        scenarioQ[s] &&
        AssociationQ[topology] &&
        Head[normalizedModel["Adjacency Matrix"]] === SparseArray
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: sparse adjacency is accepted and topology is present"
]

(* Test: user-provided sparse adjacency is preserved when Graph is also present *)
Test[
    Module[{g, model, s, normalizedModel},
        g = AdjacencyGraph[{1, 2}, {{0, 1}, {1, 0}}, DirectedEdges -> True];
        model = <|
            "Graph" -> g,
            "Vertices List" -> {1, 2},
            "Adjacency Matrix" -> SparseArray[{{1, 2} -> 3, {2, 1} -> 0}, {2, 2}],
            "Entrance Vertices and Flows" -> {{1, 10}},
            "Exit Vertices and Terminal Costs" -> {{2, 0}},
            "Switching Costs" -> <||>
        |>;
        s = makeScenario[<|"Model" -> model|>];
        normalizedModel = ScenarioData[s, "Model"];
        scenarioQ[s] &&
        Head[normalizedModel["Adjacency Matrix"]] === SparseArray &&
        normalizedModel["Adjacency Matrix"][[1, 2]] === 3
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: graph normalization preserves explicit sparse adjacency"
]

(* Test: makeScenario fails when required keys exist but topology cannot be built *)
Test[
    Module[{model, result},
        model = <|
            "Vertices List" -> {1, 2},
            "Adjacency Matrix" -> "invalid-adjacency",
            "Entrance Vertices and Flows" -> {{1, 10}},
            "Exit Vertices and Terminal Costs" -> {{2, 0}},
            "Switching Costs" -> <||>
        |>;
        result = makeScenario[<|"Model" -> model|>];
        FailureQ[result] && result["Tag"] === "ScenarioValidation"
    ]
    ,
    True
    ,
    {AdjacencyGraph::matsq, EdgeList::graph, EdgeList::graph}
    ,
    TestID -> "Scenario kernel: invalid topology shape is rejected during construction"
]

(* Test: completeScenario can be called directly on a validated scenario *)
Test[
    Module[{data, raw, validated, completed},
        data      = $testModel12;
        raw       = scenario[<|"Model" -> data|>];
        validated = validateScenario[raw];
        completed = completeScenario[validated];
        scenarioQ[completed]
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: completeScenario callable directly on validated scenario"
]

(* Test: makeScenario rejects legacy Data substitutions *)
Test[
    Module[{raw, result},
        raw = <|
            "Model" -> <|
                "Vertices List" -> {1, 2},
                "Adjacency Matrix" -> {{0, 1}, {0, 0}},
                "Entrance Vertices and Flows" -> {{1, inflowParam}},
                "Exit Vertices and Terminal Costs" -> {{2, exitCostParam}},
                "Switching Costs" -> {}
            |>,
            "Data" -> {inflowParam -> 100, exitCostParam -> 0}
        |>;
        result = makeScenario[raw];
        FailureQ[result] && result["Tag"] === "ScenarioValidation"
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: legacy Data substitutions are rejected"
]

(* Test: makeScenario fails when boundary values are symbolic *)
Test[
    Module[{raw, result},
        raw = <|
            "Model" -> <|
                "Vertices List" -> {1, 2},
                "Adjacency Matrix" -> {{0, 1}, {0, 0}},
                "Entrance Vertices and Flows" -> {{1, inflowParam}},
                "Exit Vertices and Terminal Costs" -> {{2, exitCostParam}},
                "Switching Costs" -> {}
            |>
        |>;
        result = makeScenario[raw];
        FailureQ[result] && result["Tag"] === "ScenarioValidation"
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: symbolic boundary values are rejected"
]

(* Test: makeScenario fails when switching cost values are symbolic *)
Test[
    Module[{raw, result},
        raw = <|
            "Model" -> <|
                "Vertices List" -> {1, 2},
                "Adjacency Matrix" -> {{0, 1}, {0, 0}},
                "Entrance Vertices and Flows" -> {{1, 100}},
                "Exit Vertices and Terminal Costs" -> {{2, 0}},
                "Switching Costs" -> {{1, 1, 2, c12}}
            |>
        |>;
        result = makeScenario[raw];
        FailureQ[result] && result["Tag"] === "ScenarioValidation"
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: non-numeric switching costs are rejected"
]

(* Test: Graph-only model derives Vertices List and Adjacency Matrix *)
Test[
    Module[{model, s, normalized},
        model = <|
            "Graph" -> Graph[{UndirectedEdge[1, 2], UndirectedEdge[2, 3]}],
            "Entrance Vertices and Flows" -> {{1, 10}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs" -> {}
        |>;
        s = makeScenario[<|"Model" -> model|>];
        normalized = ScenarioData[s, "Model"];
        scenarioQ[s] &&
        ListQ[normalized["Vertices List"]] &&
        MatrixQ[normalized["Adjacency Matrix"]] &&
        Length[normalized["Vertices List"]] == Length[normalized["Adjacency Matrix"]]
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: graph-only model derives topology fields"
]

(* Test: explicit Vertices List is preserved when Graph is present *)
Test[
    Module[{model, s, normalized},
        model = <|
            "Graph" -> Graph[{UndirectedEdge[1, 2], UndirectedEdge[2, 3]}],
            "Vertices List" -> {3, 2, 1},
            "Entrance Vertices and Flows" -> {{1, 10}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs" -> {}
        |>;
        s = makeScenario[<|"Model" -> model|>];
        normalized = ScenarioData[s, "Model"];
        scenarioQ[s] &&
        normalized["Vertices List"] === {3, 2, 1} &&
        MatrixQ[normalized["Adjacency Matrix"]]
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: graph normalization preserves explicit vertices order"
]

(* Test: makeScenario rejects non-integer vertex labels *)
Test[
    Module[{model, result},
        model = <|
            "Vertices List" -> {1, "2"},
            "Adjacency Matrix" -> {{0, 1}, {1, 0}},
            "Entrance Vertices and Flows" -> {{1, 10}},
            "Exit Vertices and Terminal Costs" -> {{"2", 0}},
            "Switching Costs" -> {}
        |>;
        result = makeScenario[<|"Model" -> model|>];
        FailureQ[result] && result["Tag"] === "ScenarioValidation"
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: non-integer vertices are rejected"
]

(* Test: GraphScenario with non-integer vertex labels is rejected *)
Test[
    Module[{result},
        result = GraphScenario[
            Graph[{"a" -> "b"}],
            {{"a", 1}},
            {{"b", 0}}
        ];
        FailureQ[result] && result["Tag"] === "ScenarioValidation"
    ]
    ,
    True
    ,
    TestID -> "GraphScenario-rejects-non-integer-vertices"
]

(* Test: malformed switching-cost row {1} (length != 4) must be rejected *)
Test[
    FailureQ @ makeScenario[<|
        "Model" -> <|
            "Vertices List"                   -> {1, 2},
            "Adjacency Matrix"                -> {{0, 1}, {0, 0}},
            "Entrance Vertices and Flows"     -> {{1, 1}},
            "Exit Vertices and Terminal Costs" -> {{2, 0}},
            "Switching Costs"                 -> {{1}}
        |>
    |>]
    ,
    True
    ,
    TestID -> "SwitchingCosts-malformed-row-rejected"
]

(* Test: Association-form switching costs are accepted and stored as Association *)
Test[
    Module[{s, sc},
        s = makeScenario[<|
            "Model" -> <|
                "Vertices List"                   -> {1, 2, 3, 4},
                "Adjacency Matrix"                -> {{0,1,0,0},{0,0,1,1},{0,0,0,0},{0,0,0,0}},
                "Entrance Vertices and Flows"     -> {{1, 1}},
                "Exit Vertices and Terminal Costs" -> {{3, 0}, {4, 0}},
                "Switching Costs"                 -> <|{1,2,3} -> 2, {1,2,4} -> 3|>
            |>
        |>];
        sc = ScenarioData[s, "Model"]["Switching Costs"];
        scenarioQ[s] && AssociationQ[sc] && sc[{1,2,3}] === 2
    ]
    ,
    True
    ,
    TestID -> "SwitchingCosts-association-form-accepted"
]
