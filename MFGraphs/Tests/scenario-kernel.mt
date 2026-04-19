(* Smoke tests for the scenario kernel (makeScenario, validateScenario, completeScenario, scenarioQ). *)

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
