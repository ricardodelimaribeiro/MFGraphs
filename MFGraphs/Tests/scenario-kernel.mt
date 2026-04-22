(* Smoke tests for the scenario kernel (makeScenario, validateScenario, completeScenario, scenarioQ). *)

(* Test: public symbols exist in MFGraphs` context *)
Test[
    NameQ["MFGraphs`scenario"] &&
    NameQ["MFGraphs`scenarioQ"] &&
    NameQ["MFGraphs`makeScenario"] &&
    NameQ["MFGraphs`validateScenario"] &&
    NameQ["MFGraphs`completeScenario"] &&
    NameQ["MFGraphs`ScenarioData"] &&
    NameQ["MFGraphs`ScenarioByKey"]
    ,
    True
    ,
    TestID -> "Scenario kernel: public symbols exist"
]

(* Test: makeScenario returns a typed scenario for valid input *)
Test[
    Module[{data, s},
        data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
        s = makeScenario[<|"Model" -> data|>];
        scenarioQ[s]
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: makeScenario returns typed scenario"
]

(* Test: completed scenario carries a non-empty content hash *)
Test[
    Module[{data, s, hash},
        data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
        s = makeScenario[<|"Model" -> data|>];
        hash = ScenarioData[s, "Identity"]["contentHash"];
        StringQ[hash] && StringLength[hash] > 0
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: completed scenario has contentHash"
]

(* Test: Benchmark block defaults are filled *)
Test[
    Module[{data, s, bench},
        data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
        s = makeScenario[<|"Model" -> data|>];
        bench = ScenarioData[s, "Benchmark"];
        bench["Tier"] === "core" && bench["Timeout"] === 300
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: default Benchmark tier and timeout"
]

(* Test: user-supplied Benchmark values override defaults *)
Test[
    Module[{data, s, bench},
        data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
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
        data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
        s = makeScenario[<|"Model" -> data|>];
        underlying = ScenarioData[s];
        AssociationQ[underlying] && KeyExistsQ[underlying, "Model"]
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: ScenarioData returns association"
]

(* Test: partial Identity override — user name is preserved and contentHash is always computed *)
Test[
    Module[{data, s, identity},
        data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
        s = makeScenario[<|"Model" -> data, "Identity" -> <|"name" -> "my-scenario"|>|>];
        identity = ScenarioData[s, "Identity"];
        StringQ[identity["contentHash"]] && identity["name"] === "my-scenario"
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: partial Identity override preserves name and computes hash"
]

(* Test: contentHash is stable — same input produces the same hash *)
Test[
    Module[{data, s1, s2},
        data = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
        s1 = makeScenario[<|"Model" -> data|>];
        s2 = makeScenario[<|"Model" -> data|>];
        ScenarioData[s1, "Identity"]["contentHash"] === ScenarioData[s2, "Identity"]["contentHash"]
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: contentHash is stable for identical input"
]

(* Test: completeScenario can be called directly on a validated scenario *)
Test[
    Module[{data, raw, validated, completed},
        data      = GetExampleData[12] /. {I1 -> 100, U1 -> 0};
        raw       = scenario[<|"Model" -> data|>];
        validated = validateScenario[raw];
        completed = completeScenario[validated];
        scenarioQ[completed] && StringQ[ScenarioData[completed, "Identity"]["contentHash"]]
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: completeScenario callable directly on validated scenario"
]

(* Test: makeScenario applies Data substitutions to boundary values and materializes numerics *)
Test[
    Module[{raw, s, model, entryVals, exitVals},
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
        s = MFGraphs`makeScenario[raw];
        model = MFGraphs`ScenarioData[s, "Model"];
        entryVals = Last /@ model["Entrance Vertices and Flows"];
        exitVals = Last /@ model["Exit Vertices and Terminal Costs"];
        MFGraphs`scenarioQ[s] && AllTrue[entryVals, NumericQ] && AllTrue[exitVals, NumericQ]
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: boundary values are numeric after makeScenario"
]

(* Test: makeScenario fails when boundary values remain symbolic after Data substitution *)
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
            "Data" -> {inflowParam -> 100}
        |>;
        result = MFGraphs`makeScenario[raw];
        FailureQ[result] && result["Tag"] === "ScenarioValidation"
    ]
    ,
    True
    ,
    TestID -> "Scenario kernel: non-numeric boundary values are rejected"
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

(* Test: ScenarioByKey constructs a typed scenario from input-defined boundaries *)
Test[
    Module[{s, model},
        s = ScenarioByKey[12, <|
            "Entry flows" -> <|1 -> 100|>,
            "Exit costs" -> <|4 -> 0|>,
            "Switching Costs" -> <||>
        |>];
        model = ScenarioData[s, "Model"];
        scenarioQ[s] &&
        model["Entrance Vertices and Flows"] === {{1, 100}} &&
        model["Exit Vertices and Terminal Costs"] === {{4, 0}}
    ]
    ,
    True
    ,
    TestID -> "ScenarioByKey: builds typed scenario with input-defined boundaries"
]

(* Test: ScenarioByKey ignores base example entry/exit sets *)
Test[
    Module[{s, model},
        s = ScenarioByKey[7, <|
            "Entry flows" -> <|2 -> 11|>,
            "Exit costs" -> <|4 -> 9|>,
            "Switching Costs" -> <||>
        |>];
        model = ScenarioData[s, "Model"];
        model["Entrance Vertices and Flows"] === {{2, 11}} &&
        model["Exit Vertices and Terminal Costs"] === {{4, 9}}
    ]
    ,
    True
    ,
    TestID -> "ScenarioByKey: input boundaries override base example boundaries"
]

(* Test: ScenarioByKey completes unspecified switching triples with zero *)
Test[
    Module[{s, sc},
        s = ScenarioByKey[14, <|
            "Entry flows" -> <|1 -> 10|>,
            "Exit costs" -> <|3 -> 0|>,
            "Switching Costs" -> <|{1, 2, 3} -> 5|>
        |>];
        sc = ScenarioData[s, "Model"]["Switching Costs"];
        AssociationQ[sc] &&
        KeyExistsQ[sc, {1, 2, 3}] && sc[{1, 2, 3}] == 5 &&
        KeyExistsQ[sc, {3, 2, 1}] && sc[{3, 2, 1}] == 0
    ]
    ,
    True
    ,
    TestID -> "ScenarioByKey: missing switching triples default to zero"
]

(* Test: ScenarioByKey rejects non-adjacent switching triples *)
Test[
    Module[{result},
        result = ScenarioByKey[3, <|
            "Entry flows" -> <|1 -> 1|>,
            "Exit costs" -> <|3 -> 0|>,
            "Switching Costs" -> <|{1, 3, 2} -> 7|>
        |>];
        FailureQ[result] && result["Tag"] === "ScenarioByKey"
    ]
    ,
    True
    ,
    TestID -> "ScenarioByKey: rejects non-adjacent switching triples"
]

(* Test: ScenarioByKey validates required keys *)
Test[
    Module[{result},
        result = ScenarioByKey[12, <|
            "Entry flows" -> <|1 -> 1|>,
            "Exit costs" -> <|4 -> 0|>
        |>];
        FailureQ[result] && result["Tag"] === "ScenarioByKey"
    ]
    ,
    True
    ,
    TestID -> "ScenarioByKey: missing required input keys yields Failure"
]

(* Test: ScenarioByKey reports unknown key as structured Failure *)
Test[
    Module[{result},
        result = ScenarioByKey["does-not-exist", <|
            "Entry flows" -> <|1 -> 1|>,
            "Exit costs" -> <|1 -> 0|>,
            "Switching Costs" -> <||>
        |>];
        FailureQ[result] && result["Tag"] === "ScenarioByKey"
    ]
    ,
    True
    ,
    TestID -> "ScenarioByKey: unknown key yields structured Failure"
]
