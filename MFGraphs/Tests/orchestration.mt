(* Tests for Orchestration layer *)

Test[
    NameQ["orchestrationTools`solveScenario"],
    True,
    TestID -> "Orchestration: solveScenario exists"
]

Test[
    Module[{s, result},
        s = getExampleScenario[2, {{1, 100}}, {{2, 0}}];
        result = TimeConstrained[solveScenario[s], 2, $TimedOut];
        (ListQ[result] || AssociationQ[result]) &&
        isValidSystemSolution[makeSystem[s], result]
    ],
    True,
    TestID -> "solveScenario: chains example 2 (chain 2) to a solution"
]

Test[
    Module[{s, result},
        s = getExampleScenario[2, {{1, 100}}, {{2, 0}}];
        result = TimeConstrained[solveScenario[s], 2, $TimedOut];
        (ListQ[result] || AssociationQ[result]) &&
        isValidSystemSolution[makeSystem[s], result]
    ],
    True,
    TestID -> "solveScenario: chains example 2 (chain 2) to a solution with default solver"
]

Test[
    Module[{s1, s2, results},
        s1 = getExampleScenario[2, {{1, 100}}, {{2, 0}}];
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
        s = getExampleScenario[2, {{1, 100}}, {{2, 0}}];
        result = TimeConstrained[solveScenario[s, dnfReduceSystem], 5, $TimedOut];
        (ListQ[result] || AssociationQ[result]) &&
        isValidSystemSolution[makeSystem[s], result]
    ],
    True,
    TestID -> "solveScenario: supports custom solver dnfReduceSystem"
]
