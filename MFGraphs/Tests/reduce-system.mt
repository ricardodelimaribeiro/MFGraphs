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

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2];
        sys = makeSystem[s];
        result = Quiet[reduceSystem[sys], reduceSystem::noncritical];
        FailureQ[result] && result["Tag"] === "reduceSystem"
    ],
    True,
    TestID -> "reduceSystem: non-critical congestion systems fail"
]

(* --- dnfReduceSystem tests --- *)

Test[
    NameQ["solversTools`dnfReduceSystem"],
    True,
    TestID -> "dnfReduceSystem: public symbol exists"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        (ListQ[result] && MatchQ[result, {__Rule}]) ||
        (AssociationQ[result] && KeyExistsQ[result, "Rules"])
    ],
    True,
    TestID -> "dnfReduceSystem: chain 2-exits returns rules or rules+equations"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "dnfReduceSystem: chain 2-exits solution is valid"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[8, {{1, 80.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "dnfReduceSystem: y-network solution is valid"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2];
        sys = makeSystem[s];
        result = Quiet[dnfReduceSystem[sys], dnfReduceSystem::noncritical];
        FailureQ[result] && result["Tag"] === "dnfReduceSystem"
    ],
    True,
    TestID -> "dnfReduceSystem: non-critical congestion systems fail"
]
