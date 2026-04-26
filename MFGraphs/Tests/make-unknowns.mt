(* Tests for makeUnknowns helper *)

Test[
    NameQ["unknownsTools`makeUnknowns"],
    True,
    TestID -> "makeUnknowns: public symbol exists"
]

Test[
    NameQ["unknownsTools`unknowns"] && NameQ["unknownsTools`unknownsQ"] && NameQ["unknownsTools`unknownsData"],
    True,
    TestID -> "makeUnknowns: typed unknowns symbols exist"
]

Test[
    Module[{data, s, unk},
        data = <|
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> {}
        |>;
        s = makeScenario[<|"Model" -> data|>];
        unk = makeUnknowns[s];
        unknownsQ[unk] &&
        ListQ[unknownsData[unk, "Js"]] &&
        ListQ[unknownsData[unk, "Jts"]] &&
        ListQ[unknownsData[unk, "Us"]] &&
        Length[unknownsData[unk, "Js"]] > 0 &&
        Length[unknownsData[unk, "Jts"]] > 0 &&
        Length[unknownsData[unk, "Us"]] === Length[unknownsData[unk, "Js"]] &&
        MatchQ[First[unknownsData[unk, "Js"]], _j] &&
        MatchQ[First[unknownsData[unk, "Jts"]], _j] &&
        MatchQ[First[unknownsData[unk, "Us"]], _u]
    ],
    True,
    TestID -> "makeUnknowns: builds flow, transition-flow, and value unknowns from scenario"
]

Test[
    Module[{data, s, unk, us},
        data = <|
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> {}
        |>;
        s = makeScenario[<|"Model" -> data|>];
        unk = makeUnknowns[s];
        us = unknownsData[unk, "Us"];
        !AnyTrue[
            us,
            MatchQ[#, HoldPattern[u[_, b_]] /; StringQ[b] && StringStartsQ[b, "aux"]] &
        ]
    ],
    True,
    TestID -> "makeUnknowns: boundary value unknowns use u[aux*, v] orientation only"
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
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> {}
        |>;
        unk = makeUnknowns[scenario[<|"Model" -> data|>]];
        unknownsQ[unk] &&
        Length[unknownsData[unk, "Js"]] > 0 &&
        Length[unknownsData[unk, "Jts"]] > 0 &&
        Length[unknownsData[unk, "Us"]] > 0
    ],
    True,
    TestID -> "makeUnknowns: accepts raw typed scenario head"
]

Test[
    Module[{data, customScenarioHead, unk},
        data = <|
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> {}
        |>;
        customScenarioHead = Symbol["TmpCtx`scenario"];
        unk = makeUnknowns[customScenarioHead[<|
            "Identity" -> <|"Name" -> "custom-context scenario"|>,
            "Model" -> data
        |>]];
        unknownsQ[unk] &&
        Length[unknownsData[unk, "Js"]] > 0 &&
        Length[unknownsData[unk, "Jts"]] > 0 &&
        Length[unknownsData[unk, "Us"]] > 0
    ],
    True,
    TestID -> "makeUnknowns: accepts scenario head from custom context"
]

Test[
    Module[{data, s, unk, scenarioTriples, unknownTriples},
        data = <|
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> {}
        |>;
        s = makeScenario[<|"Model" -> data|>];
        unk = makeUnknowns[s];
        scenarioTriples = scenarioData[s, "Topology"]["AuxTriples"];
        unknownTriples = unknownsData[unk, "AuxTriples"];
        ListQ[scenarioTriples] && ListQ[unknownTriples] && unknownTriples === scenarioTriples
    ],
    True,
    TestID -> "makeUnknowns: auxTriples preserve scenario topology order"
]

Test[
    Module[{data, s, unk, raw, badTriples, badUnk, sys},
        data = <|
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> {}
        |>;
        s = makeScenario[<|"Model" -> data|>];
        unk = makeUnknowns[s];
        raw = unknownsData[unk];
        If[!AssociationQ[raw], Return[False, Module]];
        badTriples = If[
            Length[raw["AuxTriples"]] > 1,
            RotateLeft[raw["AuxTriples"], 1],
            Join[raw["AuxTriples"], {{Symbol["vIn"], Symbol["vMid"], Symbol["vOut"]}}]
        ];
        badUnk = unknowns[Join[raw, <|"AuxTriples" -> badTriples|>]];
        sys = makeSystem[s, badUnk];
        FailureQ[sys] && sys["Tag"] === "makeSystem"
    ],
    True,
    TestID -> "makeSystem: mismatched auxTriples are rejected"
]

Test[
    Module[{data, s, unk, sys, scenarioTriples, systemTriples},
        data = <|
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> {}
        |>;
        s = makeScenario[<|"Model" -> data|>];
        unk = makeUnknowns[s];
        sys = makeSystem[s, unk];
        scenarioTriples = scenarioData[s, "Topology"]["AuxTriples"];
        systemTriples = systemData[sys, "AuxTriples"];
        mfgSystemQ[sys] && ListQ[systemTriples] && systemTriples === scenarioTriples
    ],
    True,
    TestID -> "makeSystem: canonical auxTriples from scenario are propagated downstream"
]

Test[
    Module[{data, s, unk, sys, entryEq, entryRules, exitValueRules, entryOutRules, exitInRules,
            B, KM, vars},
        data = <|
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> {}
        |>;
        s = makeScenario[<|"Model" -> data|>];
        unk = makeUnknowns[s];
        sys = makeSystem[s, unk];
        entryEq = systemData[sys, "EqEntryIn"];
        entryRules = systemData[sys, "RuleEntryIn"];
        exitValueRules = systemData[sys, "RuleExitValues"];
        entryOutRules = systemData[sys, "RuleEntryOut"];
        exitInRules = systemData[sys, "RuleExitFlowsIn"];
        {B, KM, vars} = getKirchhoffLinearSystem[sys];
        mfgSystemQ[sys] &&
        MemberQ[entryEq, j["auxEntry1", 1] == 10] &&
        KeyExistsQ[entryRules, j["auxEntry1", 1]] &&
        entryRules[j["auxEntry1", 1]] === 10 &&
        KeyExistsQ[exitValueRules, u["auxExit3", 3]] &&
        !KeyExistsQ[exitValueRules, u[3, "auxExit3"]] &&
        AssociationQ[entryOutRules] && entryOutRules === <||> &&
        AssociationQ[exitInRules] && exitInRules === <||> &&
        !MemberQ[vars, j["auxEntry1", 1]]
    ],
    True,
    TestID -> "makeSystem: boundary conditions include RuleEntryIn and Kirchhoff eliminates entry boundary variables"
]
