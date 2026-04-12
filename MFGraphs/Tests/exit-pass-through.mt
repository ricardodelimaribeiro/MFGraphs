(* Wolfram Language Test file *)

Test[
    Module[{data, d2e, result, asso},
        data = <|
            "Vertices List" -> {1, 2, 3},
            "Adjacency Matrix" -> {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}},
            "Entrance Vertices and Flows" -> {{1, 4}},
            "Exit Vertices and Terminal Costs" -> {{2, 10}, {3, 0}},
            "Switching Costs" -> {}
        |>;
        d2e = DataToEquations[data];
        result = CriticalCongestionSolver[d2e];
        asso = Lookup[result, "AssoCritical", <||>];
        Lookup[result, "ResultKind"] === "Success" &&
        Lookup[result, "Status"] === "Feasible" &&
        AssociationQ[asso] &&
        Lookup[asso, MFGraphs`Private`j[2, 3], Missing[]] === 4 &&
        Lookup[asso, MFGraphs`Private`j[2, ex2], Missing[]] === 0 &&
        Lookup[asso, MFGraphs`Private`j[3, ex3], Missing[]] === 4
    ],
    True,
    TestID -> "Exit pass-through: iota=4 routes all mass to the cheaper exit at v3"
]

Test[
    Module[{data, d2e, result, asso},
        data = <|
            "Vertices List" -> {1, 2, 3},
            "Adjacency Matrix" -> {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}},
            "Entrance Vertices and Flows" -> {{1, 12}},
            "Exit Vertices and Terminal Costs" -> {{2, 10}, {3, 0}},
            "Switching Costs" -> {}
        |>;
        d2e = DataToEquations[data];
        result = CriticalCongestionSolver[d2e];
        asso = Lookup[result, "AssoCritical", <||>];
        Lookup[result, "ResultKind"] === "Success" &&
        Lookup[result, "Status"] === "Feasible" &&
        AssociationQ[asso] &&
        Lookup[asso, MFGraphs`Private`j[2, 3], Missing[]] === 10 &&
        Lookup[asso, MFGraphs`Private`j[2, ex2], Missing[]] === 2 &&
        Lookup[asso, MFGraphs`Private`j[3, ex3], Missing[]] === 10
    ],
    True,
    TestID -> "Exit pass-through: iota=12 saturates v3 path and sends remainder to v2 exit"
]
