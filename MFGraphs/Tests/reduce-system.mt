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
    Module[{s, sys, result, rules, residual},
        s = getExampleScenario[7, {{1, 100.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        rules = If[ListQ[result], result, Lookup[result, "Rules", {}]];
        residual = If[AssociationQ[result], Lookup[result, "Equations", True], True];
        And[
            AssociationQ[result],
            (j[2, 4] /. rules) === 100 - j[2, 3],
            (j[3, "auxExit3"] /. rules) === j[2, 3],
            (j[4, "auxExit4"] /. rules) === 100 - j[2, 3],
            !FreeQ[residual, j[2, 3] == 55],
            isValidSystemSolution[sys, result]
        ]
    ],
    True,
    TestID -> "dnfReduceSystem: example-7 preserves edge-flow branches"
]

Test[
    Module[{s, sys, result, rules, residual},
        s = getExampleScenario[7, {{1, 100.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        rules = If[ListQ[result], result, Lookup[result, "Rules", {}]];
        residual = If[AssociationQ[result], Lookup[result, "Equations", True], True];
        And[
            AssociationQ[result],
            (u[1, 2] /. rules) === -100 + u[2, 1],
            !FreeQ[residual, u[2, 3] == 0],
            !FreeQ[residual, u[2, 4] == 10],
            isValidSystemSolution[sys, result]
        ]
    ],
    True,
    TestID -> "dnfReduceSystem: example-7 preserves value branches"
]

Test[
    Module[{s, sys, result, rules, residual},
        s = gridScenario[{3}, {{1, 5}}, {{2, 0}, {3, 10}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        rules = If[ListQ[result], result, Lookup[result, "Rules", {}]];
        residual = If[AssociationQ[result], Lookup[result, "Equations", True], True];
        And[
            AssociationQ[result],
            (j[2, 3] /. rules) === 0,
            (j[3, 2] /. rules) === 0,
            (u[2, 3] /. rules) === u[2, 3],
            TrueQ[Simplify[Implies[residual, u[2, 3] <= 10]]],
            isValidSystemSolution[sys, result]
        ]
    ],
    True,
    TestID -> "dnfReduceSystem: zero-flow edge leaves unused value constrained"
]

Test[
    Module[{s, sys, result, residual},
        s = gridScenario[{3}, {{1, 5}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        residual = If[AssociationQ[result], Lookup[result, "Equations", True], True];
        And[
            AssociationQ[result],
            !FreeQ[residual, j[2, "auxExit2"] == 5],
            !FreeQ[residual, j[2, "auxExit2"] == 0],
            isValidSystemSolution[sys, result]
        ]
    ],
    True,
    TestID -> "dnfReduceSystem: relaxed edge equation preserves chain branches"
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

Test[
    Module[{s, sys},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        sys = makeSystem[s];
        !isValidSystemSolution[
            sys,
            <|"Rules" -> {}, "Equations" -> (u[1, 2] == 0 && u[1, 2] == 1)|>
        ]
    ],
    True,
    TestID -> "isValidSystemSolution: rejects inconsistent residual"
]

(* --- booleanReduceSystem tests --- *)

Test[
    NameQ["solversTools`booleanReduceSystem"],
    True,
    TestID -> "booleanReduceSystem: public symbol exists"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = booleanReduceSystem[sys];
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "booleanReduceSystem: chain 2-exits solution is valid"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[8, {{1, 80.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = Quiet[booleanReduceSystem[sys], booleanReduceSystem::multisol];
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "booleanReduceSystem: y-network solution is valid"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2];
        sys = makeSystem[s];
        result = Quiet[booleanReduceSystem[sys], booleanReduceSystem::noncritical];
        FailureQ[result] && result["Tag"] === "booleanReduceSystem"
    ],
    True,
    TestID -> "booleanReduceSystem: non-critical congestion systems fail"
]

(* --- findInstanceSystem tests --- *)

Test[
    NameQ["solversTools`findInstanceSystem"],
    True,
    TestID -> "findInstanceSystem: public symbol exists"
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
        result = findInstanceSystem[sys];
        MatchQ[result, {__Rule}] && isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "findInstanceSystem: chain 1-exit returns valid rules"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = findInstanceSystem[sys];
        MatchQ[result, {__Rule}] && isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "findInstanceSystem: chain 2-exits solution is valid"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[8, {{1, 80.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = findInstanceSystem[sys];
        MatchQ[result, {__Rule}] && isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "findInstanceSystem: y-network solution is valid"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2];
        sys = makeSystem[s];
        result = Quiet[findInstanceSystem[sys], findInstanceSystem::noncritical];
        FailureQ[result] && result["Tag"] === "findInstanceSystem"
    ],
    True,
    TestID -> "findInstanceSystem: non-critical congestion systems fail"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[8, {{1, 80.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = findInstanceSystem[sys, "Timeout" -> 0];
        AssociationQ[result] &&
        KeyExistsQ[result, "Rules"] &&
        Lookup[result, "Equations", Missing["KeyAbsent", "Equations"]] === False
    ],
    True,
    TestID -> "findInstanceSystem: timeout returns rules plus false residual"
]
