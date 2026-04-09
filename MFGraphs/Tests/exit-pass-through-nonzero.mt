(* Wolfram Language Test file *)

Test[
    Module[{vertices, adjacency, graph, auxiliaryGraph, auxVertices, auxTriples,
            switchingCosts, data, d2e, result, asso, j23, j2ex2, j3ex3},
        vertices = {1, 2, 3};
        adjacency = {{0, 1, 0}, {1, 0, 1}, {0, 1, 0}};
        graph = AdjacencyGraph[vertices, adjacency, DirectedEdges -> False];
        auxiliaryGraph = EdgeAdd[
            graph,
            {
                UndirectedEdge[en1, 1],
                UndirectedEdge[2, ex2],
                UndirectedEdge[3, ex3]
            }
        ];
        auxVertices = VertexList[auxiliaryGraph];
        auxTriples = Flatten[
            Insert[#, 2] & /@ Permutations[AdjacencyList[auxiliaryGraph, #], {2}] & /@ auxVertices,
            1
        ];
        (* Force the nonzero-switching preprocessing path while keeping a consistent metric. *)
        switchingCosts = Append[#, 1] & /@ auxTriples;

        data = <|
            "Vertices List" -> vertices,
            "Adjacency Matrix" -> adjacency,
            "Entrance Vertices and Flows" -> {{1, 12}},
            "Exit Vertices and Terminal Costs" -> {{2, 10}, {3, 0}},
            "Switching Costs" -> switchingCosts
        |>;
        d2e = DataToEquations[data];
        result = CriticalCongestionSolver[d2e];
        asso = Lookup[result, "AssoCritical", <||>];
        j23 = Lookup[asso, MFGraphs`Private`j[2, 3], Missing[]];
        j2ex2 = Lookup[asso, MFGraphs`Private`j[2, ex2], Missing[]];
        j3ex3 = Lookup[asso, MFGraphs`Private`j[3, ex3], Missing[]];

        Lookup[result, "ResultKind"] === "Success" &&
        Lookup[result, "Status"] === "Feasible" &&
        AssociationQ[asso] &&
        NumericQ[j23] && NumericQ[j2ex2] && NumericQ[j3ex3] &&
        j23 > 0 &&
        Chop[j23 + j2ex2 - 12] === 0 &&
        Chop[j3ex3 - j23] === 0
    ],
    True,
    TestID -> "Exit pass-through: nonzero switching costs keep model feasible and preserve flow through v3 path"
]
