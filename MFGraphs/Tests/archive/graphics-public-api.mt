(* Minimal regression tests for the public graphics API extracted in Issue #25. *)

Test[
    NameQ["MFGraphs`NetworkGraphPlot"] &&
    NameQ["MFGraphs`SolutionFlowPlot"] &&
    NameQ["MFGraphs`ExitFlowPlot"]
    ,
    True
    ,
    TestID -> "Graphics API: public symbols exist"
]

Test[
    Module[{d2e, networkPlot, flowPlot, exitPlot},
        d2e = <|
            "graph" -> Graph[{DirectedEdge[1, 2]}],
            "entryVertices" -> {1},
            "exitVertices" -> {2},
            "edgeList" -> {DirectedEdge[1, 2]},
            "SignedFlows" -> <|{1, 2} -> 5.|>,
            "RuleBalanceGatheringFlows" -> <||>,
            "RuleExitFlowsIn" -> <||>,
            "RuleEntryOut" -> <||>
        |>;
        networkPlot = NetworkGraphPlot[d2e, "Network smoke test"];
        flowPlot = SolutionFlowPlot[d2e, <||>, "Flow smoke test"];
        exitPlot = ExitFlowPlot[<|2 -> 5.|>, "Exit smoke test"];
        MatchQ[
            {Head[networkPlot], Head[flowPlot], Head[exitPlot]},
            {Graph, Graph, Graphics}
        ]
    ]
    ,
    True
    ,
    TestID -> "Graphics API: plotting helpers are callable"
]
