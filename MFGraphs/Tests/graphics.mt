(* Tests for Graphics helpers *)

Test[
    NameQ["graphics`scenarioTopologyPlot"] && NameQ["graphics`mfgSolutionPlot"],
    True,
    TestID -> "Graphics: public symbols exist"
]

Test[
    Module[{s, sys, sol, topoPlot, solPlot},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
        topoPlot = scenarioTopologyPlot[s, sys];
        solPlot = mfgSolutionPlot[s, sys, sol];
        MatchQ[topoPlot, _Graph] && MatchQ[solPlot, _Graph]
    ],
    True,
    TestID -> "Graphics: topology and solution plots return Graph"
]

Test[
    Module[{s, sys, sol, edges},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
        edges = EdgeList[mfgSolutionPlot[s, sys, sol]];
        MemberQ[edges, DirectedEdge[1, 2]]
    ],
    True,
    TestID -> "Graphics: positive j[1,2] displays as DirectedEdge[1,2]"
]
