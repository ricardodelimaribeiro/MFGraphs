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
        Module[{res},
            d2e = DataToEquations[data];
            res = Quiet[MonotoneSolverFromData[data, opts]];
            {res, Lookup[res, "ComparableFlowVector", $Failed]}
        ]
    ];

(* Shape regression: a simple path solve should return the standardized association. *)
Test[
    Module[{data, res, flows},
        data = GetExampleData[3] /. {I1 -> 80, U1 -> 0};
        {res, flows} = monotoneSolveWithFlows[data, "ResidualTolerance" -> 10^-6, "MaxTime" -> 10, "MaxSteps" -> 2000];
        AssociationQ[res] && Lookup[res, "ResultKind"] === "Success" && AssociationQ[res["AssoMonotone"]]
    ]
    ,
    True
    ,
    TestID -> "Monotone solver: path case returns the standardized result association"
]

(* Flow regression: the 2-edge path should recover 80 units of net flow on each
   interior edge. The 10^-4 tolerance allows for NDSolve/integration noise. *)
Test[
    Module[{data, res, flows},
        data = GetExampleData[3] /. {I1 -> 80, U1 -> 0};
        {res, flows} = monotoneSolveWithFlows[data, "ResidualTolerance" -> 10^-6, "MaxTime" -> 10, "MaxSteps" -> 2000];
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
        {res1, flows1} = monotoneSolveWithFlows[data, "UseCachedProjection" -> True, "ResidualTolerance" -> 10^-6, "MaxTime" -> 10, "MaxSteps" -> 2000];
        {res2, flows2} = monotoneSolveWithFlows[data, "UseCachedProjection" -> False, "ResidualTolerance" -> 10^-6, "MaxTime" -> 10, "MaxSteps" -> 2000];
        AssociationQ[res1] && AssociationQ[res2] &&
        Lookup[res1, "ResultKind"] === "Success" &&
        Lookup[res2, "ResultKind"] === "Success" &&
        Max[Abs[flows1 - flows2]] < 10^-5
    ]
    ,
    True
    ,
    TestID -> "Monotone solver: cached and uncached path flows agree"
]

(* Critical-case regression: on the asymmetric Y-network the Monotone solver
   should now use the fast quadratic critical path and recover the same net
   edge flows as CriticalCongestionSolver. *)
Test[
    Module[{data, d2e, crit, res, flows},
        data = GetExampleData[7] /. {I1 -> 80, U1 -> 0, U2 -> 10};
        d2e = DataToEquations[data];
        crit = Quiet[CriticalCongestionSolver[d2e]];
        {res, flows} = monotoneSolveWithFlows[data, "ResidualTolerance" -> 10^-6, "MaxTime" -> 20, "MaxSteps" -> 5000];
        Lookup[res, "ResultKind"] === "Success" &&
        Lookup[res, "Feasibility"] === "Feasible" &&
        Lookup[res["Convergence"], "StopReason"] === "QuadraticCriticalSolve" &&
        Max[Abs[flows - crit["ComparableFlowVector"]]] < 10^-6
    ]
    ,
    True
    ,
    TestID -> "Monotone solver: asymmetric Y-network matches the critical congestion solution"
]

(* Switching regression: on the 2-exit switching network the quadratic critical
   path should recover the same split as CriticalCongestionSolver. *)
Test[
    Module[{data, d2e, crit, res, flows},
        data = GetExampleData[8] /. {I1 -> 100, U1 -> 0, U2 -> 10, S1 -> 2, S2 -> 3, S3 -> 2, S4 -> 1, S5 -> 3, S6 -> 1};
        d2e = DataToEquations[data];
        crit = Quiet[CriticalCongestionSolver[d2e]];
        {res, flows} = monotoneSolveWithFlows[data, "ResidualTolerance" -> 10^-6, "MaxTime" -> 20, "MaxSteps" -> 5000];
        Lookup[res, "ResultKind"] === "Success" &&
        Lookup[res, "Feasibility"] === "Feasible" &&
        Lookup[res["Convergence"], "StopReason"] === "QuadraticCriticalSolve" &&
        Max[Abs[flows - crit["ComparableFlowVector"]]] < 10^-6
    ]
    ,
    True
    ,
    TestID -> "Monotone solver: switching Y-network matches the critical congestion solution"
]

(* Paper example regression: keep the published Scenario 1 data wired into the
   repo without making the fast suite solve the full paper-scale problem. *)
Test[
    Module[{data, d2e},
        data = GetExampleData["HRF Scenario 1"] /. {I1 -> 100, I2 -> 100, U1 -> 0, U2 -> 0};
        d2e = DataToEquations[data];
        Lookup[data, "Entrance Vertices and Flows"] === {{1, 100}, {9, 100}} &&
        Lookup[data, "Exit Vertices and Terminal Costs"] === {{8, 0}, {10, 0}} &&
        Length[Lookup[d2e, "edgeList", {}]] === 15
    ]
    ,
    True
    ,
    TestID -> "Monotone solver: HRF Scenario 1 data stays wired as a 15-edge paper example"
]

(* Fast-path regression: the quadratic critical solve should ignore tiny ODE
   budgets because it does not use the iterative ODE path at all. *)
Test[
    Module[{data, res, flows},
        data = GetExampleData[7] /. {I1 -> 80, U1 -> 0, U2 -> 10};
        {res, flows} = monotoneSolveWithFlows[data, "ResidualTolerance" -> 10^-12, "MaxTime" -> 20, "MaxSteps" -> 1];
        Lookup[res, "ResultKind"] === "Success" &&
        Lookup[res["Convergence"], "StopReason"] === "QuadraticCriticalSolve"
    ]
    ,
    True
    ,
    TestID -> "Monotone solver: quadratic critical path bypasses the ODE step budget"
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
        MatchQ[
            Quiet[MonotoneSolverFromData[data]],
            _Association?(Lookup[#, "ResultKind", None] === "Degenerate" && Lookup[#, "Message", None] === "DegenerateCase" &)
        ]
    ]
    ,
    True
    ,
    TestID -> "Monotone solver: explicit empty-flow input is degenerate"
]
