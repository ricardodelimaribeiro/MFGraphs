(* Tests for makeSymbolicUnknowns helper *)

Test[
    NameQ["unknownsTools`makeSymbolicUnknowns"],
    True,
    TestID -> "makeSymbolicUnknowns: public symbol exists"
]

Test[
    NameQ["unknownsTools`symbolicUnknowns"] && NameQ["unknownsTools`symbolicUnknownsQ"] && NameQ["unknownsTools`symbolicUnknownsData"],
    True,
    TestID -> "makeSymbolicUnknowns: typed symbolicUnknowns symbols exist"
]

Test[
    NameQ["unknownsTools`unknown"] && NameQ["unknownsTools`unknowns"] &&
        !NameQ["unknownsTools`makeUnknowns"] &&
        !NameQ["unknownsTools`unknownsQ"] &&
        !NameQ["unknownsTools`unknownsData"],
    True,
    TestID -> "symbolicUnknowns: legacy helper names are not public aliases"
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
        unk = makeSymbolicUnknowns[s];
        symbolicUnknownsQ[unk] &&
        ListQ[symbolicUnknownsData[unk, "Js"]] &&
        ListQ[symbolicUnknownsData[unk, "Jts"]] &&
        ListQ[symbolicUnknownsData[unk, "Us"]] &&
        Length[symbolicUnknownsData[unk, "Js"]] > 0 &&
        Length[symbolicUnknownsData[unk, "Jts"]] > 0 &&
        Length[symbolicUnknownsData[unk, "Us"]] === Length[symbolicUnknownsData[unk, "Js"]] &&
        MatchQ[First[symbolicUnknownsData[unk, "Js"]], _j] &&
        MatchQ[First[symbolicUnknownsData[unk, "Jts"]], _j] &&
        MatchQ[First[symbolicUnknownsData[unk, "Us"]], _u]
    ],
    True,
    TestID -> "makeSymbolicUnknowns: builds flow, transition-flow, and value variables from scenario"
]

Test[
    Module[{data, s, unk, js, us, jts},
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
        unk = makeSymbolicUnknowns[s];
        js = symbolicUnknownsData[unk, "Js"];
        us = symbolicUnknownsData[unk, "Us"];
        jts = symbolicUnknownsData[unk, "Jts"];
        MemberQ[js, j["auxEntry1", 1]] &&
        MemberQ[js, j[1, 2]] &&
        MemberQ[us, u["auxEntry1", 1]] &&
        MemberQ[us, u[1, 2]] &&
        MemberQ[jts, j["auxEntry1", 1, 2]] &&
        MemberQ[jts, j[1, 2, 3]]
    ],
    True,
    TestID -> "makeSymbolicUnknowns: preserves exact j and u variable forms"
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
        unk = makeSymbolicUnknowns[s];
        us = symbolicUnknownsData[unk, "Us"];
        !AnyTrue[
            us,
            MatchQ[#, HoldPattern[u[_, b_]] /; StringQ[b] && StringStartsQ[b, "aux"]] &
        ]
    ],
    True,
    TestID -> "makeSymbolicUnknowns: boundary value variables use u[aux*, v] orientation only"
]

Test[
    Module[{result},
        result = makeSymbolicUnknowns[<|"Model" -> <||>|>];
        FailureQ[result]
    ],
    True,
    TestID -> "makeSymbolicUnknowns: rejects non-scenario input"
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
        unk = makeSymbolicUnknowns[scenario[<|"Model" -> data|>]];
        symbolicUnknownsQ[unk] &&
        Length[symbolicUnknownsData[unk, "Js"]] > 0 &&
        Length[symbolicUnknownsData[unk, "Jts"]] > 0 &&
        Length[symbolicUnknownsData[unk, "Us"]] > 0
    ],
    True,
    TestID -> "makeSymbolicUnknowns: accepts raw typed scenario head"
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
        unk = makeSymbolicUnknowns[customScenarioHead[<|
            "Identity" -> <|"Name" -> "custom-context scenario"|>,
            "Model" -> data
        |>]];
        symbolicUnknownsQ[unk] &&
        Length[symbolicUnknownsData[unk, "Js"]] > 0 &&
        Length[symbolicUnknownsData[unk, "Jts"]] > 0 &&
        Length[symbolicUnknownsData[unk, "Us"]] > 0
    ],
    True,
    TestID -> "makeSymbolicUnknowns: accepts scenario head from custom context"
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
        unk = makeSymbolicUnknowns[s];
        scenarioTriples = scenarioData[s, "Topology"]["AuxTriples"];
        unknownTriples = symbolicUnknownsData[unk, "AuxTriples"];
        ListQ[scenarioTriples] && ListQ[unknownTriples] && unknownTriples === scenarioTriples
    ],
    True,
    TestID -> "makeSymbolicUnknowns: auxTriples preserve scenario topology order"
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
        unk = makeSymbolicUnknowns[s];
        raw = symbolicUnknownsData[unk];
        If[!AssociationQ[raw], Return[False, Module]];
        badTriples = If[
            Length[raw["AuxTriples"]] > 1,
            RotateLeft[raw["AuxTriples"], 1],
            Join[raw["AuxTriples"], {{Symbol["vIn"], Symbol["vMid"], Symbol["vOut"]}}]
        ];
        badUnk = symbolicUnknowns[Join[raw, <|"AuxTriples" -> badTriples|>]];
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
        unk = makeSymbolicUnknowns[s];
        sys = makeSystem[s, unk];
        scenarioTriples = scenarioData[s, "Topology"]["AuxTriples"];
        systemTriples = systemData[sys, "AuxTriples"];
        mfgSystemQ[sys] && ListQ[systemTriples] && systemTriples === scenarioTriples
    ],
    True,
    TestID -> "makeSystem: canonical auxTriples from scenario are propagated downstream"
]

Test[
    Module[{data, s, sys},
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
        sys = makeSystem[s];
        mfgSystemQ[sys] &&
        MemberQ[systemData[sys, "Js"], j["auxEntry1", 1]] &&
        MemberQ[systemData[sys, "Us"], u["auxEntry1", 1]]
    ],
    True,
    TestID -> "makeSystem: automatically builds symbolic unknowns"
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
        unk = makeSymbolicUnknowns[s];
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
