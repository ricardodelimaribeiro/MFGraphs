(* Tests for makeUnknowns helper *)

Test[
    NameQ["MFGraphs`makeUnknowns"],
    True,
    TestID -> "makeUnknowns: public symbol exists"
]

Test[
    NameQ["MFGraphs`unknowns"] && NameQ["MFGraphs`unknownsQ"] && NameQ["MFGraphs`UnknownsData"],
    True,
    TestID -> "makeUnknowns: typed unknowns symbols exist"
]

Test[
    Module[{data, s, unk},
        data = <|
            "Vertices List" -> {1, 2, 3},
            "Adjacency Matrix" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entrance Vertices and Flows" -> {{1, 10}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs" -> {}
        |>;
        s = makeScenario[<|"Model" -> data|>];
        unk = makeUnknowns[s];
        unknownsQ[unk] &&
        ListQ[UnknownsData[unk, "js"]] &&
        ListQ[UnknownsData[unk, "jts"]] &&
        ListQ[UnknownsData[unk, "us"]] &&
        Length[UnknownsData[unk, "js"]] > 0 &&
        Length[UnknownsData[unk, "jts"]] > 0 &&
        Length[UnknownsData[unk, "us"]] === Length[UnknownsData[unk, "js"]] &&
        MatchQ[First[UnknownsData[unk, "js"]], _j] &&
        MatchQ[First[UnknownsData[unk, "jts"]], _j] &&
        MatchQ[First[UnknownsData[unk, "us"]], _u]
    ],
    True,
    TestID -> "makeUnknowns: builds flow, transition-flow, and value unknowns from scenario"
]

Test[
    Module[{result},
        result = makeUnknowns[<|"Model" -> <||>|>];
        FailureQ[result]
    ],
    True,
    TestID -> "makeUnknowns: rejects non-scenario input"
]

Test[
    Module[{data, unk},
        data = <|
            "Vertices List" -> {1, 2, 3},
            "Adjacency Matrix" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entrance Vertices and Flows" -> {{1, 10}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs" -> {}
        |>;
        unk = makeUnknowns[scenario[<|"Model" -> data|>]];
        unknownsQ[unk] &&
        Length[UnknownsData[unk, "js"]] > 0 &&
        Length[UnknownsData[unk, "jts"]] > 0 &&
        Length[UnknownsData[unk, "us"]] > 0
    ],
    True,
    TestID -> "makeUnknowns: accepts raw typed scenario head"
]

Test[
    Quiet[
        Module[{data, gscenario, unk},
            data = <|
                "Vertices List" -> {1, 2, 3},
                "Adjacency Matrix" -> {
                    {0, 1, 0},
                    {1, 0, 1},
                    {0, 1, 0}
                },
                "Entrance Vertices and Flows" -> {{1, 10}},
                "Exit Vertices and Terminal Costs" -> {{3, 0}},
                "Switching Costs" -> {}
            |>;
            gscenario = Symbol["Global`scenario"];
            unk = makeUnknowns[gscenario[<|
                "Identity" -> <|"Name" -> "shadowed scenario"|>,
                "Model" -> data
            |>]];
            unknownsQ[unk] &&
            Length[UnknownsData[unk, "js"]] > 0 &&
            Length[UnknownsData[unk, "jts"]] > 0 &&
            Length[UnknownsData[unk, "us"]] > 0
        ],
        {Global`scenario::shdw}
    ],
    True,
    TestID -> "makeUnknowns: accepts scenario head from shadowed context"
]

Test[
    Module[{data, s, unk, scenarioTriples, unknownTriples},
        data = <|
            "Vertices List" -> {1, 2, 3},
            "Adjacency Matrix" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entrance Vertices and Flows" -> {{1, 10}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs" -> {}
        |>;
        s = makeScenario[<|"Model" -> data|>];
        unk = makeUnknowns[s];
        scenarioTriples = ScenarioData[s, "Topology"]["AuxTriples"];
        unknownTriples = UnknownsData[unk, "auxTriples"];
        ListQ[scenarioTriples] && ListQ[unknownTriples] && unknownTriples === scenarioTriples
    ],
    True,
    TestID -> "makeUnknowns: auxTriples preserve scenario topology order"
]

Test[
    Module[{data, s, unk, raw, badTriples, badUnk, sys},
        data = <|
            "Vertices List" -> {1, 2, 3},
            "Adjacency Matrix" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entrance Vertices and Flows" -> {{1, 10}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs" -> {}
        |>;
        s = makeScenario[<|"Model" -> data|>];
        unk = makeUnknowns[s];
        raw = UnknownsData[unk];
        If[!AssociationQ[raw], Return[False, Module]];
        badTriples = If[
            Length[raw["auxTriples"]] > 1,
            RotateLeft[raw["auxTriples"], 1],
            Join[raw["auxTriples"], {{Symbol["vIn"], Symbol["vMid"], Symbol["vOut"]}}]
        ];
        badUnk = unknowns[Join[raw, <|"auxTriples" -> badTriples|>]];
        sys = makeSystem[s, badUnk];
        FailureQ[sys] && sys["Tag"] === "makeSystem"
    ],
    True,
    TestID -> "makeSystem: mismatched auxTriples are rejected"
]

Test[
    Module[{data, s, unk, sys, scenarioTriples, systemTriples},
        data = <|
            "Vertices List" -> {1, 2, 3},
            "Adjacency Matrix" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entrance Vertices and Flows" -> {{1, 10}},
            "Exit Vertices and Terminal Costs" -> {{3, 0}},
            "Switching Costs" -> {}
        |>;
        s = makeScenario[<|"Model" -> data|>];
        unk = makeUnknowns[s];
        sys = makeSystem[s, unk];
        scenarioTriples = ScenarioData[s, "Topology"]["AuxTriples"];
        systemTriples = SystemData[sys, "auxTriples"];
        mfgSystemQ[sys] && ListQ[systemTriples] && systemTriples === scenarioTriples
    ],
    True,
    TestID -> "makeSystem: canonical auxTriples from scenario are propagated downstream"
]
