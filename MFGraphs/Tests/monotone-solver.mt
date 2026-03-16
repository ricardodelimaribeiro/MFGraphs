(* Wolfram Language Test file *)
(* Solver-facing regression tests for MonotoneSolver *)

(* Local helper for solver-facing regressions.
   V[x, edge] is the edge potential term in the Hamiltonian H.
   alpha[edge] is the congestion exponent used in H and Cost.
   g[m, edge] is the density interaction term in H.
   We fix all three with Block so Cost and the induced monotone field are
   reproducible even if the surrounding session overrides these symbols.
   The helper returns net interior edge flows j[a,b] - j[b,a] rather than the
   full j[...] association because auxiliary and transition variables are a
   noisier regression target. A negative signed value means net flow in the
   reverse orientation, not negative mass. *)
monotoneSolveWithFlows[data_, opts___] :=
    Block[{V = Function[{x, edge}, 0],
           alpha = Function[edge, 1],
           g = Function[{m, edge}, -1/m^2]},
        Module[{d2e, res, flowExprs},
            d2e = DataToEquations[data];
            res = Quiet[MonotoneSolverFromData[data, opts]];
            flowExprs = (d2e["SignedFlows"] /@ (List @@@ d2e["edgeList"])) /.
                Join[d2e["RuleBalanceGatheringFlows"], d2e["RuleExitFlowsIn"], d2e["RuleEntryOut"]];
            {res, If[AssociationQ[res], N[flowExprs /. res], $Failed]}
        ]
    ];

(* Shape regression: a simple path solve should return an association keyed by j[...]. *)
Test[
    Module[{data, res, flows},
        data = GetExampleData[3] /. {I1 -> 80, U1 -> 0};
        {res, flows} = monotoneSolveWithFlows[data, "TimeSteps" -> 20];
        AssociationQ[res]
    ]
    ,
    True
    ,
    TestID -> "Monotone solver: path case returns an association"
]

(* Flow regression: the 2-edge path should recover 80 units of net flow on each
   interior edge. The 10^-4 tolerance allows for NDSolve/integration noise. *)
Test[
    Module[{data, res, flows},
        data = GetExampleData[3] /. {I1 -> 80, U1 -> 0};
        {res, flows} = monotoneSolveWithFlows[data, "TimeSteps" -> 20];
        Max[Abs[flows - {80., 80.}]] < 10^-4
    ]
    ,
    True
    ,
    TestID -> "Monotone solver: path case gives 80 flow on both edges"
]

(* Cache regression: cached and uncached gradient projection should agree on the
   same path problem. The 10^-5 tolerance is tighter because both solves target
   the same net interior flows. *)
Test[
    Module[{data, res1, res2, flows1, flows2},
        data = GetExampleData[3] /. {I1 -> 80, U1 -> 0};
        {res1, flows1} = monotoneSolveWithFlows[data, "UseCachedGradient" -> True, "TimeSteps" -> 20];
        {res2, flows2} = monotoneSolveWithFlows[data, "UseCachedGradient" -> False, "TimeSteps" -> 20];
        AssociationQ[res1] && AssociationQ[res2] && Max[Abs[flows1 - flows2]] < 10^-5
    ]
    ,
    True
    ,
    TestID -> "Monotone solver: cached and uncached path flows agree"
]

(* Degenerate-input regression: explicit raw data with no flow variables should
   short-circuit before ODE solving. This uses raw data because no built-in
   example exercises the no-flow path directly. *)
Test[
    Module[{data},
        data = <|
            "Vertices List" -> {1},
            "Adjacency Matrix" -> {{0}},
            "Entrance Vertices and Flows" -> {},
            "Exit Vertices and Terminal Costs" -> {},
            "Switching Costs" -> {}
        |>;
        Quiet[MonotoneSolverFromData[data]]
    ]
    ,
    <|"Message" -> "Degenerate case"|>
    ,
    TestID -> "Monotone solver: explicit empty-flow input is degenerate"
]
