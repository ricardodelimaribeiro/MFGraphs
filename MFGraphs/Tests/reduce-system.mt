(* Tests for ReduceSystem *)

Test[
    NameQ["MFGraphs`ReduceSystem"],
    True,
    TestID -> "ReduceSystem: public symbol exists"
]

Test[
    Module[{data, s, sys, result},
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
        sys = makeSystem[s];
        result = ReduceSystem[sys];
        result =!= False
    ],
    True,
    TestID -> "ReduceSystem: chain 1-exit yields non-False solution"
]

Test[
    Module[{s, sys, result},
        s = GridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = ReduceSystem[sys];
        result =!= False
    ],
    True,
    TestID -> "ReduceSystem: chain 2-exits yields non-False solution"
]
