(* Tests for Orchestration layer *)

Test[
    NameQ["orchestrationTools`solveScenario"],
    True,
    TestID -> "Orchestration: solveScenario exists"
]

Test[
    Module[{s, result},
        s = gridScenario[{2}, {{1, 100}}, {{2, 0}}];
        result = TimeConstrained[solveScenario[s], 2, $TimedOut];
        (ListQ[result] || AssociationQ[result]) &&
        isValidSystemSolution[makeSystem[s], result]
    ],
    True,
    TestID -> "solveScenario: chains example 2 (chain 2) to a solution"
]

Test[
    Module[{s, result},
        s = gridScenario[{2}, {{1, 100}}, {{2, 0}}];
        result = TimeConstrained[solveScenario[s], 2, $TimedOut];
        (ListQ[result] || AssociationQ[result]) &&
        isValidSystemSolution[makeSystem[s], result]
    ],
    True,
    TestID -> "solveScenario: chains example 2 (chain 2) to a solution with default solver"
]

Test[
    Module[{s1, s2, results},
        s1 = gridScenario[{2}, {{1, 100}}, {{2, 0}}];
        s2 = getExampleScenario[3, {{1, 100}}, {{3, 0}}]; (* Example 3 is another simple chain *)
        results = TimeConstrained[solveScenario[{s1, s2}], 5, $TimedOut];
        ListQ[results] && Length[results] === 2 &&
        isValidSystemSolution[makeSystem[s1], results[[1]]] &&
        isValidSystemSolution[makeSystem[s2], results[[2]]]
    ],
    True,
    TestID -> "solveScenario: handles multi-population (list of scenarios)"
]

Test[
    Module[{s, result},
        s = gridScenario[{2}, {{1, 100}}, {{2, 0}}];
        result = TimeConstrained[solveScenario[s, dnfReduceSystem], 5, $TimedOut];
        (ListQ[result] || AssociationQ[result]) &&
        isValidSystemSolution[makeSystem[s], result]
    ],
    True,
    TestID -> "solveScenario: supports custom solver dnfReduceSystem"
]

Test[
    NameQ["orchestrationTools`SolveMFG"],
    True,
    TestID -> "Orchestration: SolveMFG exists"
]

Test[
    Module[{s, solveScenarioResult, solveMFGResult},
        s = getExampleScenario[3, {{1, 10}}, {{3, 0}}];
        solveScenarioResult = TimeConstrained[solveScenario[s], 5, $TimedOut];
        solveMFGResult = TimeConstrained[SolveMFG[s], 5, $TimedOut];
        (ListQ[solveMFGResult] || AssociationQ[solveMFGResult]) &&
        solutionResultKind[solveMFGResult] === solutionResultKind[solveScenarioResult]
    ],
    True,
    TestID -> "SolveMFG: accepts typed scenario and produces same result kind as solveScenario"
]

Test[
    Module[{rawAssoc, s, solveScenarioResult, solveMFGResult},
        rawAssoc = <|"Model" -> <|
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {{0,1,0},{0,0,1},{0,0,0}},
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> {}
        |>|>;
        s = makeScenario[rawAssoc];
        solveScenarioResult = TimeConstrained[solveScenario[s], 5, $TimedOut];
        solveMFGResult = TimeConstrained[SolveMFG[rawAssoc], 5, $TimedOut];
        (ListQ[solveMFGResult] || AssociationQ[solveMFGResult]) &&
        solutionResultKind[solveMFGResult] === solutionResultKind[solveScenarioResult]
    ],
    True,
    TestID -> "SolveMFG: legacy association input still works after scenario dispatch added"
]
