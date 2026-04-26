(* Tests for Graphics helpers *)

Test[
    NameQ["MFGraphs`ScenarioTopologyPlot"] && NameQ["MFGraphs`MFGSolutionPlot"],
    True,
    TestID -> "Graphics: public symbols exist"
]

Test[
    Module[{s, sys, sol, topoPlot, solPlot},
        s = GridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = ReduceSystem[sys];
        topoPlot = ScenarioTopologyPlot[s, sys];
        solPlot = MFGSolutionPlot[s, sys, sol];
        MatchQ[topoPlot, _Graph] && MatchQ[solPlot, _Graph]
    ],
    True,
    TestID -> "Graphics: topology and solution plots return Graph"
]

Test[
    Module[{s, sys, sol, edges},
        s = GridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = ReduceSystem[sys];
        edges = EdgeList[MFGSolutionPlot[s, sys, sol]];
        MemberQ[edges, DirectedEdge[1, 2]]
    ],
    True,
    TestID -> "Graphics: positive j[1,2] displays as DirectedEdge[1,2]"
]
