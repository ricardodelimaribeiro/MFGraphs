(* Wolfram Language package *)
(* Iterative solver for the non-critical (general) congestion case *)

Clear[H, Cost];

(* --- Options --- *)

Options[NonLinearSolver] = {"MaxIterations" -> 15, "Tolerance" -> 0};

(* --- NonLinearSolver: main iterative solver --- *)

NonLinearSolver::usage =
    "NonLinearSolver[Eqs] takes an association resulting from DataToEquations and returns an approximation to the solution of the non-critical congestion case with alpha = value, specified by alpha[edge_] := value.
Options: \"MaxIterations\" (default 15), \"Tolerance\" (default 0). When Tolerance > 0, iteration stops early when the infinity-norm change in flow variables between consecutive steps falls below the given tolerance.";

NonLinearSolver[Eqs_, OptionsPattern[]] :=
    Module[ {AssoCritical, PreEqs = Eqs, AssoNonCritical, NonCriticalList, js,
             MaxIter = OptionValue["MaxIterations"], tol = OptionValue["Tolerance"]},
        If[ KeyExistsQ[PreEqs, "AssoCritical"],
            (* If there is already an approximation for the non-congestion case, use it *)
            AssoNonCritical = Lookup[PreEqs, "AssoNonCritical", PreEqs["AssoCritical"]],
            PreEqs = MFGPreprocessing[PreEqs];
            AssoNonCritical = AssociationThread[Lookup[PreEqs, "js", {}], 0 Lookup[PreEqs, "js", {}]]
        ];
        js = Lookup[PreEqs, "js", {}];
        NonCriticalList = If[ tol > 0 && js =!= {},
            (* Tolerance-based stopping: compare consecutive flow vectors *)
            NestWhileList[
                nonLinearStep[PreEqs],
                AssoNonCritical,
                (Norm[N[Values[KeyTake[#2, js]] - Values[KeyTake[#1, js]]], Infinity] > tol)&,
                2, MaxIter
            ],
            (* Default: exact fixed-point check (stops when consecutive results are identical) *)
            FixedPointList[nonLinearStep[PreEqs], AssoNonCritical, MaxIter]
        ];
        MFGPrint["Iterated ", Length[NonCriticalList]-1, " times out of ", MaxIter];
        AssoCritical = Lookup[PreEqs, "AssoCritical", NonCriticalList[[2]]];
        AssoNonCritical = NonCriticalList // Last;
        Join[PreEqs, Association[{"AssoCritical" -> AssoCritical, "AssoNonCritical" -> AssoNonCritical}]]
    ];

(* --- nonLinearStep: single iteration --- *)

nonLinearStep[Eqs_][approxSol_] :=
    Module[ {approxJs, approx, js, Nrhs, Nlhs, Newlhs, Newrhs},
        js = Lookup[Eqs, "js", $Failed];
        Nrhs = Lookup[Eqs, "Nrhs", $Failed];
        Nlhs = Lookup[Eqs, "Nlhs", $Failed];
        approxJs = KeyTake[approxSol, js];
        approx = MFGSystemSolver[Eqs][approxJs];
        Newlhs = N[Nlhs/.approx];
        Newrhs = Nrhs/.approx;
        MFGPrint[Newlhs,"\n",Newrhs];
        MFGPrint[Style["Max error for non-linear solution: ", Bold, Blue], Norm[Newlhs-Newrhs, Infinity]];
        approx
    ];

(* --- IsNonLinearSolution: validation --- *)

IsNonLinearSolution::usage =
"IsNonLinearSolution[Eqs] extracts AssoNonCritical and checks equations and inequalities. The right and left hand sides of the nonlinear equations are shown with the sup-norm of the difference."

IsNonLinearSolution[Eqs_] :=
    Reap@Module[ {EqEntryIn, EqValueAuxiliaryEdges, IneqSwitchingByVertex, AltOptCond,
        EqBalanceSplittingFlows, AltFlows, AltTransitionFlows,
        IneqJs, IneqJts, Nrhs, Nlhs,
        styleBlue, styleRed, styleGreen, assoc, allSatisfied},
        assoc = Lookup[Eqs, "AssoNonCritical", Print["Needs \"AssoNonCritical\" solution."]; Sow[<||>]];
        styleBlue = Style[#, Bold, Blue]&;
        styleRed = Style[#, Bold, Red]&;
        styleGreen = Style[#, Bold, Darker@Green]&;
        EqEntryIn = And @@ Lookup[Eqs, "EqEntryIn", $Failed];
        AltOptCond = Lookup[Eqs, "AltOptCond", $Failed];
        AltFlows = Lookup[Eqs, "AltFlows", $Failed];
        AltTransitionFlows = Lookup[Eqs, "AltTransitionFlows", $Failed];
        IneqJs = Lookup[Eqs, "IneqJs", $Failed];
        IneqJts = Lookup[Eqs, "IneqJts", $Failed];
        IneqSwitchingByVertex = And @@ Lookup[Eqs, "IneqSwitchingByVertex", $Failed];
        EqBalanceSplittingFlows = Lookup[Eqs, "EqBalanceSplittingFlows", $Failed];
        EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", $Failed];
        Nrhs = Lookup[Eqs, "Nrhs", $Failed];
        Nlhs = Lookup[Eqs, "Nlhs", $Failed];
        allSatisfied = (EqEntryIn&&EqValueAuxiliaryEdges&&AltOptCond&&
            EqBalanceSplittingFlows&&AltFlows&&AltTransitionFlows&&
            IneqJs&&IneqJts&&IneqSwitchingByVertex)/.assoc;
        If[ allSatisfied,
            Print["All restrictions are ", styleGreen@"True"],
            Print["At least one of the restrictions is ", styleRed@allSatisfied];
            Print[styleBlue@"EqEntryIn: ", EqEntryIn/.assoc, "\n", EqEntryIn];
            Print[styleBlue@"EqValueAuxiliaryEdges: ", EqValueAuxiliaryEdges/.assoc, "\n", EqValueAuxiliaryEdges];
            Print[styleBlue@"AltOptCond: ", AltOptCond/.assoc, "\n", AltOptCond];
            Print[styleBlue@"EqBalanceSplittingFlows: ", EqBalanceSplittingFlows/.assoc, "\n", EqBalanceSplittingFlows];
            Print[styleBlue@"AltFlows: ", AltFlows/.assoc, "\n",AltFlows];
            Print[styleBlue@"AltTransitionFlows: ", AltTransitionFlows/.assoc, "\n", AltTransitionFlows];
            Print[styleBlue@"IneqJs: ", IneqJs/.assoc, "\n", IneqJs];
            Print[styleBlue@"IneqJts: ", IneqJts/.assoc, "\n", IneqJts];
            Print[styleBlue@"IneqSwitchingByVertex: ", IneqSwitchingByVertex/.assoc, "\n", IneqSwitchingByVertex]
        ];
        
        Print[styleBlue@"Nlhs: ", Nlhs/.assoc, "\n", Nlhs];
        Print[styleBlue@"Nrhs: ", Nrhs/.assoc, "\n", Nrhs];
        Print[styleBlue@"Max error for non-linear solution: ", Norm[N[Nlhs/.assoc]-(Nrhs/.assoc), Infinity]];
        N@assoc
    ];

(* --- Hamiltonian and related functions --- *)

(* Default congestion exponent *)
alpha[edge_] := 1

(* Default interaction potential *)
g[m_, edge_] := -1/m^2

(* Potential function - users should override V[x, edge] for specific problems *)
(* Default: V = Function[{x, edge}, W[x, A]] where W[y,a] = a Sin[2 Pi (y+1/4)]^2 *)

H::usage =
"H[xi, p, m, edge] is the Hamiltonian function for the edges.
edge is a directed edge from the graph."
H = Function[{xi, p, m, edge}, p^2/(2 m^alpha[edge]) + V[xi, edge] - g[m, edge]];

U::usage =
"U[x, edge, Eqs, sol] computes the value function at position x on the given edge."
U[x_?NumericQ , edge_, Eqs_Association, sol_] :=
    Module[ {jay, uT},
        jay = Eqs["SignedFlows"][List@@edge]/.sol;
        uT = u@@Reverse@(List@@edge);
        uT = uT/.sol;
        If[ PossibleZeroQ[jay],
            uT,
            uT - jay NIntegrate[M[jay, y, edge]^(alpha[edge] - 1), {y, 0, x}]
        ]
    ]

M[j_?NumericQ, x_?NumericQ, edge_] :=
    If[ PossibleZeroQ[j],
        0.,
        Values@First@FindRoot[H[x, -j m^(alpha[edge] - 1), m, edge], {m, 1}]
    ];

IntegratedMass[j_?NumericQ, edge_] :=
    If[ PossibleZeroQ[j],
        0,
        j NIntegrate[M[j, x, edge]^(alpha[edge] - 1), {x, 0, 1}] // Quiet
    ];

Cost[j_, edge_] :=
    IntegratedMass[j, edge];

(* --- Interpolation-based M for faster Cost evaluation --- *)

PrecomputeM::usage =
"PrecomputeM[jMin, jMax, edge, nPoints] precomputes M[j,x,edge] on a grid and returns
an InterpolatingFunction. This avoids per-point FindRoot calls during NIntegrate.
nPoints controls the grid resolution (default 50).";

PrecomputeM[jMin_?NumericQ, jMax_?NumericQ, edge_, nPoints_Integer:50] :=
  Module[{jVals, xVals, mGrid, flatData},
    jVals = Subdivide[jMin, jMax, nPoints];
    xVals = Subdivide[0., 1., nPoints];
    mGrid = Table[
      If[PossibleZeroQ[jv], 0.,
        Quiet @ Check[
          Values @ First @ FindRoot[H[xv, -jv m^(alpha[edge]-1), m, edge], {m, 1}],
          0.
        ]
      ],
      {jv, jVals}, {xv, xVals}
    ];
    flatData = Flatten[
      Table[
        {{jVals[[i]], xVals[[k]]}, mGrid[[i, k]]},
        {i, Length@jVals}, {k, Length@xVals}
      ], 1
    ];
    Interpolation[flatData, InterpolationOrder -> 3]
  ];

FastIntegratedMass::usage =
"FastIntegratedMass[interpM, j, edge] computes IntegratedMass using a precomputed interpolation of M.
interpM should be the result of PrecomputeM.";

FastIntegratedMass[interpM_, j_?NumericQ, edge_] :=
  If[PossibleZeroQ[j], 0,
    j NIntegrate[interpM[j, x]^(alpha[edge]-1), {x, 0, 1}] // Quiet
  ];

(* --- Plotting utilities --- *)

PlotMassDensity[Eqs_Association, string_String, pair_List] :=
    Module[ {jays, sol, edge},
        sol = Eqs[string];
        jays = Eqs["SignedFlows"][pair] /. sol;
        edge = UndirectedEdge@@pair;
        Plot[M[jays /. sol, x, edge], {x, 0, 1},
         PlotLabel -> edge,
         GridLines -> Automatic]
    ];

PlotMassDensities[Eqs_, string_] :=
    PlotMassDensity[Eqs, string, #] & /@ Lookup[Eqs, "pairs", {}]

PlotValueFunction[Eqs_Association, string_String, pair_List] :=
    Module[ {sol, edge},
        sol = Eqs[string];
        edge = UndirectedEdge@@pair;
        Plot[U[x, edge, Eqs, sol], {x, 0, 1}, PlotLabel -> edge,
         GridLines -> Automatic]
    ];

PlotValueFunctions[Eqs_, string_] :=
    PlotValueFunction[Eqs, string, #] & /@ Lookup[Eqs, "pairs", {}]
