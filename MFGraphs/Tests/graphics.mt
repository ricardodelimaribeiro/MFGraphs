(* Tests for Graphics helpers *)

Test[
    NameQ["graphics`scenarioTopologyPlot"] && NameQ["graphics`mfgSolutionPlot"],
    True,
    TestID -> "Graphics: public symbols exist"
]

Test[
    Module[{names, globalNames},
        names = {"scenarioTopologyPlot", "mfgSolutionPlot"};
        globalNames = ("Global`" <> #) & /@ names;
        Scan[If[NameQ[#], Remove[#]] &, globalNames];
        Needs["MFGraphs`"];
        And @@ (NameQ /@ names)
    ],
    True,
    TestID -> "Graphics: unqualified symbols are available after Needs"
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

Test[
    Module[{s, sys, sol, edges},
        s = gridScenario[{3}, {{3, 120}}, {{1, 0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
        edges = EdgeList[mfgSolutionPlot[s, sys, sol]];
        MemberQ[edges, DirectedEdge[2, 1]]
    ],
    True,
    TestID -> "Graphics: negative net flow flips display direction to DirectedEdge[2,1]"
]

Test[
    Module[{s, sys, topoPlot, solPlot},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        topoPlot = scenarioTopologyPlot[s, sys, "Topology custom title"];
        solPlot = mfgSolutionPlot[s, sys, <||>, "Solution custom title"];
        !FreeQ[AbsoluteOptions[topoPlot, PlotLabel], "Topology custom title"] &&
        !FreeQ[AbsoluteOptions[solPlot, PlotLabel], "Solution custom title"]
    ],
    True,
    TestID -> "Graphics: custom titles propagate to PlotLabel"
]

Test[
    Module[{s, sys, sol, rules},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        sol = reduceSystem[sys];
        rules = Lookup[sol, "Rules", {}];
        MatchQ[mfgSolutionPlot[s, sys, rules], _Graph]
    ],
    True,
    TestID -> "Graphics: mfgSolutionPlot accepts raw rule list solution"
]

Test[
    Module[{s, sys, graph},
        s = gridScenario[{3}, {{1, 120}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        graph = mfgSolutionPlot[s, sys, <||>];
        MatchQ[graph, _Graph] && Length[EdgeList[graph]] > 0
    ],
    True,
    TestID -> "Graphics: mfgSolutionPlot handles missing Rules association without failure"
]
