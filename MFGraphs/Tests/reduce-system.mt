(* Tests for reduceSystem *)

Test[
    NameQ["solversTools`reduceSystem"],
    True,
    TestID -> "reduceSystem: public symbol exists"
]

Test[
    Module[{data, s, sys, result},
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
        result = reduceSystem[sys];
        result =!= False
    ],
    True,
    TestID -> "reduceSystem: chain 1-exit yields non-False solution"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = reduceSystem[sys];
        result =!= False
    ],
    True,
    TestID -> "reduceSystem: chain 2-exits yields non-False solution"
]

Test[
    Module[{s, sys, entryVals, exitVals},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        entryVals = Values @ Normal @ systemData[sys, "RuleEntryIn"];
        exitVals = Values @ Normal @ systemData[sys, "RuleExitValues"];
        FreeQ[Join[entryVals, exitVals], _Real]
    ],
    True,
    TestID -> "reduceSystem: boundary rules are exactified (no Real coefficients)"
]

Test[
    Module[{s, sys, exitRules},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        exitRules = systemData[sys, "RuleExitValues"];
        KeyExistsQ[exitRules, u["auxExit2", 2]] &&
        KeyExistsQ[exitRules, u["auxExit3", 3]] &&
        !KeyExistsQ[exitRules, u[2, "auxExit2"]] &&
        !KeyExistsQ[exitRules, u[3, "auxExit3"]]
    ],
    True,
    TestID -> "reduceSystem: RuleExitValues use u[auxExit, vertex] orientation"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 10.0}, {3, 0.0}}];
        sys = makeSystem[s];
        result = reduceSystem[sys];
        result =!= False
    ],
    True,
    TestID -> "reduceSystem: chain with costs {2,10},{3,0} remains solvable"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = reduceSystem[sys];
        result =!= False
    ],
    True,
    TestID -> "reduceSystem: chain with costs {2,0},{3,10} remains solvable via inequality+complementarity"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 10.5}}, {{3, 0.25}}];
        sys = makeSystem[s];
        result = reduceSystem[sys];
        result =!= False
    ],
    True,
    TestID -> "reduceSystem: non-integer decimal boundaries solve after exactification"
]
