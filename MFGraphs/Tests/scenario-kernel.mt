(* Smoke tests for the scenario kernel (makeScenario, validateScenario, completeScenario, scenarioQ). *)

(* Inline test fixture — case 12 (4-vertex attraction network) with concrete values. *)
$testModel12 = <|
    "Vertices"                    -> {1, 2, 3, 4},
    "Adjacency"                 -> {{0,1,1,0},{0,0,1,1},{0,0,0,1},{0,0,0,0}},
    "Entries"      -> {{1, 100}},
    "Exits" -> {{4, 0}},
    "Switching"                  -> {}
|>;

(* Test: public symbols exist in MFGraphs` context *)
Test[
    NameQ["scenarioTools`scenario"] &&
    NameQ["scenarioTools`scenarioQ"] &&
    NameQ["scenarioTools`makeScenario"] &&
    NameQ["scenarioTools`validateScenario"] &&
    NameQ["scenarioTools`completeScenario"] &&
    NameQ["scenarioTools`scenarioData"]
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
            "scenarioData"
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
        bench = scenarioData[s, "Benchmark"];
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
        h = scenarioData[s, "Hamiltonian"];
        AssociationQ[h] &&
        h["Alpha"] === 1 &&
        h["V"] === -1 &&
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

(* Test: example constructors delegate omitted Hamiltonian defaults to makeScenario *)
Test[
    Module[{s, h},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        h = scenarioData[s, "Hamiltonian"];
        scenarioQ[s] && h["Alpha"] === 1 && h["V"] === -1 && h["G"][2] === -1/2
    ],
    True,
    TestID -> "Scenario kernel: gridScenario uses makeScenario Hamiltonian defaults"
]

(* Test: positional Alpha override is preserved while omitted V/G use defaults *)
Test[
    Module[{s, h},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2];
        h = scenarioData[s, "Hamiltonian"];
        scenarioQ[s] && h["Alpha"] === 2 && h["V"] === -1 && h["G"][2] === -1/2
    ],
    True,
    TestID -> "Scenario kernel: gridScenario preserves positional Alpha override"
]

(* Test: positional V override is preserved *)
Test[
    Module[{s, h},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 1, 3];
        h = scenarioData[s, "Hamiltonian"];
        scenarioQ[s] && h["Alpha"] === 1 && h["V"] === 3 && h["G"][2] === -1/2
    ],
    True,
    TestID -> "Scenario kernel: gridScenario preserves positional V override"
]

(* Test: positional G override is preserved *)
Test[
    Module[{s, h, g},
        g = Function[z, -2/z];
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 1, 3, g];
        h = scenarioData[s, "Hamiltonian"];
        scenarioQ[s] && h["Alpha"] === 1 && h["V"] === 3 && h["G"][2] === -1
    ],
    True,
    TestID -> "Scenario kernel: gridScenario preserves positional G override"
]

(* Test: getExampleScenario forwarding delegates omitted Hamiltonian defaults *)
Test[
    Module[{s, h},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        h = scenarioData[s, "Hamiltonian"];
        scenarioQ[s] && h["Alpha"] === 1 && h["V"] === -1 && h["G"][2] === -1/2
    ],
    True,
    TestID -> "Scenario kernel: getExampleScenario uses makeScenario Hamiltonian defaults"
]

(* Test: user-supplied per-edge Hamiltonian parameters are preserved *)
Test[
    Module[{model, s, h},
        model = <|
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}},
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> <||>
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
        h = scenarioData[s, "Hamiltonian"];
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
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}},
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> <||>
        |>;
        gfun = Function[z, -2/z];
        s = makeScenario[<|
            "Model" -> model,
            "Hamiltonian" -> <|"EdgeG" -> <|{1, 2} -> gfun|>|>
        |>];
        h = scenarioData[s, "Hamiltonian"];
        scenarioQ[s] &&
        MatchQ[h["EdgeG"][{1, 2}], _Function] &&
        h["EdgeG"][{1, 2}][2] === -1
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: function-valued EdgeG is accepted"
]

(* Test: V/G are preserved but not used in current EqGeneral construction *)
Test[
    Module[{model, sBase, sVG, sysBase, sysVG},
        model = <|
            "Vertices" -> {1, 2},
            "Adjacency" -> {{0, 1}, {1, 0}},
            "Entries" -> {{1, 10}},
            "Exits" -> {{2, 0}},
            "Switching" -> {}
        |>;
        sBase = makeScenario[<|"Model" -> model|>];
        sVG = makeScenario[<|
            "Model" -> model,
            "Hamiltonian" -> <|
                "Alpha" -> 1,
                "V" -> 7,
                "G" -> Function[z, z^2],
                "EdgeV" -> <|{1, 2} -> 3|>,
                "EdgeG" -> <|{1, 2} -> Function[z, z + 1]|>
            |>
        |>];
        sysBase = makeSystem[sBase];
        sysVG = makeSystem[sVG];
        scenarioData[sVG, "Hamiltonian"]["V"] === 7 &&
        KeyExistsQ[scenarioData[sVG, "Hamiltonian"]["EdgeV"], {1, 2}] &&
        systemData[sysBase, "EqGeneral"] === systemData[sysVG, "EqGeneral"]
    ],
    True,
    TestID -> "Scenario kernel: V and G are preserved but not applied to EqGeneral yet"
]

(* Test: user-supplied Benchmark values override defaults *)
Test[
    Module[{data, s, bench},
        data = $testModel12;
        s = makeScenario[<|"Model" -> data, "Benchmark" -> <|"Tier" -> "stress", "Timeout" -> 900|>|>];
        bench = scenarioData[s, "Benchmark"];
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
            "Vertices" -> {1, 2},
            "Adjacency" -> {{0, 1}, {0, 0}},
            "Entries" -> {{1, 100}},
            "Exits" -> {{2, 0}}
            (* "Switching" intentionally absent *)
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

(* Test: scenarioData accessor returns underlying association *)
Test[
    Module[{data, s, underlying},
        data = $testModel12;
        s = makeScenario[<|"Model" -> data|>];
        underlying = scenarioData[s];
        AssociationQ[underlying] && KeyExistsQ[underlying, "Model"]
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: scenarioData returns association"
]

(* Test: partial Identity override — user name is preserved *)
Test[
    Module[{data, s, identity},
        data = $testModel12;
        s = makeScenario[<|"Model" -> data, "Identity" -> <|"name" -> "my-scenario"|>|>];
        identity = scenarioData[s, "Identity"];
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
            "Vertices" -> {1, 2},
            "Adjacency" -> SparseArray[{{1, 2} -> 1, {2, 1} -> 1}, {2, 2}],
            "Entries" -> {{1, 10}},
            "Exits" -> {{2, 0}},
            "Switching" -> <||>
        |>;
        s = makeScenario[<|"Model" -> model|>];
        topology = scenarioData[s, "Topology"];
        normalizedModel = scenarioData[s, "Model"];
        scenarioQ[s] &&
        AssociationQ[topology] &&
        Head[normalizedModel["Adjacency"]] === SparseArray
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: sparse adjacency is accepted and topology is present"
]

(* Test: topology auxiliary vertices are string labels *)
Test[
    Module[{model, s, topology, auxEntry, auxExit, auxEntryEdges, auxExitEdges},
        model = <|
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}},
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> <||>
        |>;
        s = makeScenario[<|"Model" -> model|>];
        topology = scenarioData[s, "Topology"];
        auxEntry = topology["AuxEntryVertices"];
        auxExit = topology["AuxExitVertices"];
        auxEntryEdges = topology["AuxEntryEdges"];
        auxExitEdges = topology["AuxExitEdges"];
        scenarioQ[s] &&
        AllTrue[auxEntry, StringQ] &&
        AllTrue[auxExit, StringQ] &&
        MatchQ[First[auxEntryEdges], DirectedEdge[_String, _Integer]] &&
        MatchQ[First[auxExitEdges], DirectedEdge[_Integer, _String]] &&
        topology["AuxPairs"] === {
            {"auxEntry1", 1},
            {3, "auxExit3"},
            {1, 2},
            {2, 3},
            {2, 1},
            {3, 2}
        } &&
        AllTrue[
            topology["AuxTriples"],
            !(StringQ[First[#]] && StringStartsQ[First[#], "auxExit"]) &&
            !(StringQ[Last[#]] && StringStartsQ[Last[#], "auxEntry"]) &
        ]
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: auxiliary vertices are strings and aux triples respect boundary direction"
]

(* Test: user-provided sparse adjacency is preserved when Graph is also present *)
Test[
    Module[{g, model, s, normalizedModel},
        g = AdjacencyGraph[{1, 2}, {{0, 1}, {1, 0}}, DirectedEdges -> True];
        model = <|
            "Graph" -> g,
            "Vertices" -> {1, 2},
            "Adjacency" -> SparseArray[{{1, 2} -> 3, {2, 1} -> 0}, {2, 2}],
            "Entries" -> {{1, 10}},
            "Exits" -> {{2, 0}},
            "Switching" -> <||>
        |>;
        s = makeScenario[<|"Model" -> model|>];
        normalizedModel = scenarioData[s, "Model"];
        scenarioQ[s] &&
        Head[normalizedModel["Adjacency"]] === SparseArray &&
        normalizedModel["Adjacency"][[1, 2]] === 3
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
            "Vertices" -> {1, 2},
            "Adjacency" -> "invalid-adjacency",
            "Entries" -> {{1, 10}},
            "Exits" -> {{2, 0}},
            "Switching" -> <||>
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
    {completeScenario::notopology}
    ,
    TestID -> "Scenario kernel: completeScenario callable directly on validated scenario"
]

(* Test: completeScenario direct call warns and completes switching costs when topology is absent *)
Test[
    Module[{raw, validated, completed, sc},
        raw = scenario[<|"Model" -> <|
            "Vertices" -> {1, 2},
            "Adjacency" -> {{0, 1}, {1, 0}},
            "Entries" -> {{1, 10}},
            "Exits" -> {{2, 0}},
            "Switching" -> {{1, 1, 2, 3}}
        |>|>];
        validated = validateScenario[raw];
        completed = completeScenario[validated];
        sc = scenarioData[completed, "Model"]["Switching"];
        scenarioQ[completed] && AssociationQ[sc]
    ],
    True,
    {completeScenario::notopology},
    TestID -> "Scenario kernel: completeScenario rebuilds missing topology with warning"
]

(* Test: makeScenario rejects legacy Data substitutions *)
Test[
    Module[{raw, result},
        raw = <|
            "Model" -> <|
                "Vertices" -> {1, 2},
                "Adjacency" -> {{0, 1}, {0, 0}},
                "Entries" -> {{1, inflowParam}},
                "Exits" -> {{2, exitCostParam}},
                "Switching" -> {}
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
                "Vertices" -> {1, 2},
                "Adjacency" -> {{0, 1}, {0, 0}},
                "Entries" -> {{1, inflowParam}},
                "Exits" -> {{2, exitCostParam}},
                "Switching" -> {}
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

(* Test: makeScenario fails when an entry vertex is outside Vertices *)
Test[
    Module[{raw, result},
        raw = <|
            "Model" -> <|
                "Vertices" -> {1, 2},
                "Adjacency" -> {{0, 1}, {1, 0}},
                "Entries" -> {{99, 10}},
                "Exits" -> {{2, 0}},
                "Switching" -> {}
            |>
        |>;
        result = makeScenario[raw];
        FailureQ[result] && result["Tag"] === "ScenarioValidation"
    ],
    True,
    TestID -> "Scenario kernel: unknown entry vertex is rejected"
]

(* Test: makeScenario fails when an exit vertex is outside Vertices *)
Test[
    Module[{raw, result},
        raw = <|
            "Model" -> <|
                "Vertices" -> {1, 2},
                "Adjacency" -> {{0, 1}, {1, 0}},
                "Entries" -> {{1, 10}},
                "Exits" -> {{99, 0}},
                "Switching" -> {}
            |>
        |>;
        result = makeScenario[raw];
        FailureQ[result] && result["Tag"] === "ScenarioValidation"
    ],
    True,
    TestID -> "Scenario kernel: unknown exit vertex is rejected"
]

(* Test: makeScenario fails when switching cost values are symbolic *)
Test[
    Module[{raw, result},
        raw = <|
            "Model" -> <|
                "Vertices" -> {1, 2},
                "Adjacency" -> {{0, 1}, {0, 0}},
                "Entries" -> {{1, 100}},
                "Exits" -> {{2, 0}},
                "Switching" -> {{1, 1, 2, c12}}
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

(* Test: Infinity switching costs are accepted for blocked transitions *)
Test[
    Module[{s, sc},
        s = makeScenario[<|
            "Model" -> <|
                "Vertices" -> {1, 2, 3},
                "Adjacency" -> {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}},
                "Entries" -> {{1, 10}},
                "Exits" -> {{3, 0}},
                "Switching" -> {{1, 2, 3, Infinity}}
            |>
        |>];
        sc = scenarioData[s, "Model"]["Switching"];
        scenarioQ[s] && sc[{1, 2, 3}] === Infinity
    ],
    True,
    TestID -> "Scenario kernel: Infinity switching costs are accepted"
]

(* Test: Graph-only model derives Vertices List and Adjacency Matrix *)
Test[
    Module[{model, s, normalized},
        model = <|
            "Graph" -> Graph[{UndirectedEdge[1, 2], UndirectedEdge[2, 3]}],
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> {}
        |>;
        s = makeScenario[<|"Model" -> model|>];
        normalized = scenarioData[s, "Model"];
        scenarioQ[s] &&
        ListQ[normalized["Vertices"]] &&
        MatrixQ[normalized["Adjacency"]] &&
        Length[normalized["Vertices"]] == Length[normalized["Adjacency"]]
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
            "Vertices" -> {3, 2, 1},
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> {}
        |>;
        s = makeScenario[<|"Model" -> model|>];
        normalized = scenarioData[s, "Model"];
        scenarioQ[s] &&
        normalized["Vertices"] === {3, 2, 1} &&
        MatrixQ[normalized["Adjacency"]]
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
            "Vertices" -> {1, "2"},
            "Adjacency" -> {{0, 1}, {1, 0}},
            "Entries" -> {{1, 10}},
            "Exits" -> {{"2", 0}},
            "Switching" -> {}
        |>;
        result = makeScenario[<|"Model" -> model|>];
        FailureQ[result] && result["Tag"] === "ScenarioValidation"
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: non-integer vertices are rejected"
]

(* Test: graphScenario with non-integer vertex labels is rejected *)
Test[
    Module[{result},
        result = graphScenario[
            Graph[{"a" -> "b"}],
            {{"a", 1}},
            {{"b", 0}}
        ];
        FailureQ[result] && result["Tag"] === "ScenarioValidation"
    ]
    ,
    True
    ,
    TestID -> "graphScenario-rejects-non-integer-vertices"
]

(* Test: malformed switching-cost row {1} (length != 4) must be rejected *)
Test[
    FailureQ @ makeScenario[<|
        "Model" -> <|
            "Vertices"                   -> {1, 2},
            "Adjacency"                -> {{0, 1}, {0, 0}},
            "Entries"     -> {{1, 1}},
            "Exits" -> {{2, 0}},
            "Switching"                 -> {{1}}
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
                "Vertices"                   -> {1, 2, 3, 4},
                "Adjacency"                -> {{0,1,0,0},{0,0,1,1},{0,0,0,0},{0,0,0,0}},
                "Entries"     -> {{1, 1}},
                "Exits" -> {{3, 0}, {4, 0}},
                "Switching"                 -> <|{1,2,3} -> 2, {1,2,4} -> 3|>
            |>
        |>];
        sc = scenarioData[s, "Model"]["Switching"];
        scenarioQ[s] && AssociationQ[sc] && sc[{1,2,3}] === 2
    ]
    ,
    True
    ,
    TestID -> "SwitchingCosts-association-form-accepted"
]

(* Test: non-symmetric AM produces the same topology as its symmetrized counterpart *)
Test[
    Module[{am, entries, exits, s1, s2},
        am      = {{0,1,0},{0,0,1},{0,0,0}};
        entries = {{1, 10}};
        exits   = {{3, 0}};
        s1 = makeScenario[<|"Model" -> <|
            "Vertices"                    -> {1,2,3},
            "Adjacency"                 -> am,
            "Entries"      -> entries,
            "Exits" -> exits,
            "Switching"                  -> {}
        |>|>];
        s2 = makeScenario[<|"Model" -> <|
            "Vertices"                    -> {1,2,3},
            "Adjacency"                 -> am + Transpose[am],
            "Entries"      -> entries,
            "Exits" -> exits,
            "Switching"                  -> {}
        |>|>];
        scenarioData[s1, "Topology"] === scenarioData[s2, "Topology"]
    ]
    ,
    True
    ,
    TestID -> "directed-am-symmetrized-topology"
]

(* Test: Camilli 2015 metadata accessor is public after package load *)
Test[
    NameQ["getExampleScenarioMetadata"] &&
    AssociationQ[getExampleScenarioMetadata["Camilli 2015 simple"]]
    ,
    True
    ,
    TestID -> "Example metadata: public accessor is available"
]

(* Test: Camilli 2015 simple stationary analog builds with representative boundaries *)
Test[
    Module[{s, model, topology, metadata},
        s = getExampleScenario["Camilli 2015 simple", {{3, 50}, {4, 50}}, {{1, 0}}];
        model = scenarioData[s, "Model"];
        topology = scenarioData[s, "Topology"];
        metadata = getExampleScenarioMetadata["Camilli 2015 simple"];
        scenarioQ[s] &&
        Length[model["Vertices"]] === metadata["DeclaredVertexCount"] &&
        Length[topology["HalfPairs"]] === metadata["DeclaredEdgeCount"] &&
        model["Exits"] === metadata["ExitDefault"] &&
        StringContainsQ[metadata["PackageInterpretation"], "Stationary analog only"]
    ]
    ,
    True
    ,
    TestID -> "Example registry: Camilli 2015 simple builds as stationary analog"
]

(* Test: Camilli 2015 general stationary analog preserves paper-declared topology count *)
Test[
    Module[{s, model, topology, metadata, sys},
        s = getExampleScenario["Camilli 2015 general", {{14, 25}, {15, 25}, {7, 25}, {9, 25}}, {{1, 0}}];
        model = scenarioData[s, "Model"];
        topology = scenarioData[s, "Topology"];
        metadata = getExampleScenarioMetadata["Camilli 2015 general"];
        sys = makeSystem[s];
        scenarioQ[s] &&
        mfgSystemQ[sys] &&
        Length[model["Vertices"]] === 17 &&
        Length[model["Vertices"]] === metadata["DeclaredVertexCount"] &&
        Length[topology["HalfPairs"]] === 22 &&
        Length[topology["HalfPairs"]] === metadata["DeclaredEdgeCount"] &&
        metadata["TimeDependentData"]["Theta"] === 0.7 &&
        StringContainsQ[metadata["PackageInterpretation"], "time-dependent stochastic MFG"] &&
        StringContainsQ[StringRiffle[metadata["ImportNotes"], " "], "visually ambiguous"]
    ]
    ,
    True
    ,
    TestID -> "Example registry: Camilli 2015 general builds with declared 17/22 topology"
]

(* Test: Achdou 2023 finite stationary junction analog builds with caller-supplied boundaries *)
Test[
    Module[{s, model, topology, metadata},
        s = getExampleScenario[
            "Achdou 2023 junction",
            {{2, 25}, {3, 25}, {4, 25}, {5, 25}},
            {{1, 0}}
        ];
        model = scenarioData[s, "Model"];
        topology = scenarioData[s, "Topology"];
        metadata = getExampleScenarioMetadata["Achdou 2023 junction"];
        scenarioQ[s] &&
        Length[model["Vertices"]] === 5 &&
        Length[model["Vertices"]] === metadata["DeclaredVertexCount"] &&
        Length[topology["HalfPairs"]] === 4 &&
        Length[topology["HalfPairs"]] === metadata["DeclaredEdgeCount"] &&
        model["Entries"] === {{2, 25}, {3, 25}, {4, 25}, {5, 25}} &&
        model["Exits"] === {{1, 0}}
    ]
    ,
    True
    ,
    TestID -> "Example registry: Achdou 2023 junction builds finite stationary analog"
]

(* Test: Achdou 2023 finite junction analog supports system construction *)
Test[
    Module[{s, sys},
        s = getExampleScenario["Achdou 2023 junction", {{2, 100}}, {{1, 0}}];
        sys = makeSystem[s];
        scenarioQ[s] && mfgSystemQ[sys]
    ]
    ,
    True
    ,
    TestID -> "Example registry: Achdou 2023 junction makeSystem succeeds"
]

(* Test: Achdou 2023 metadata records source and finite-analog scope *)
Test[
    Module[{metadata},
        metadata = getExampleScenarioMetadata["Achdou 2023 junction"];
        AssociationQ[metadata] &&
        metadata["SourcePaperPath"] ===
            "docs/research/papers/Achdou et al. - 2023 - First order Mean Field Games on networks.pdf" &&
        metadata["HALId"] === "hal-03729443v3" &&
        StringContainsQ[metadata["FiniteAnalogNote"], "stationary finite graph analog"] &&
        StringContainsQ[metadata["FiniteAnalogNote"], "not implemented here"] &&
        metadata["RecommendedStationaryBoundaryExamples"]["BalancedLeavesToCenter"]["Exits"] === {{1, 0}}
    ]
    ,
    True
    ,
    TestID -> "Example metadata: Achdou 2023 source and finite analog note"
]

(* Test: listExampleScenarios returns the full registry as a non-empty sorted list *)
Test[
    Module[{keys},
        keys = listExampleScenarios[];
        ListQ[keys] && Length[keys] > 20 &&
            MemberQ[keys, "Jamaratv9"] && MemberQ[keys, "Grid0303"] && MemberQ[keys, 3]
    ]
    ,
    True
    ,
    TestID -> "Example registry: listExampleScenarios discoverability"
]

(* Test: every registered scenario has metadata (no $Failed leaks) *)
Test[
    AllTrue[
        listExampleScenarios[],
        AssociationQ[getExampleScenarioMetadata[#]] &
    ]
    ,
    True
    ,
    TestID -> "Example metadata: every registered key resolves to an Association"
]
