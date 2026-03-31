(* Wolfram Language package *)
(*
   NonLinearSolver: Iterative solver for the non-critical (general) congestion case.
   
   Theoretical Basis:
   - "First-order mean-field games on networks and Wardrop equilibrium" (Al Saleh et al., 2024)
   
   This solver handles Mean Field Games where the congestion strength alpha is not necessarily 
   equal to 1. It employs a fixed-point iteration (or NestWhile with tolerance) to resolve 
   non-linear coupling between density and value functions.
*)

(* --- Public API declarations --- *)

NonLinearSolver::usage =
    "NonLinearSolver[Eqs] takes an association resulting from DataToEquations and returns an approximation to the solution of the non-critical congestion case with alpha = value, specified by alpha[edge_] := value.
Options: \"MaxIterations\" (default 15), \"Tolerance\" (default 0), \"ReturnShape\" \
(default \"Legacy\"; use \"Standard\" to add normalized solver-result keys), \
\"PotentialFunction\", \"CongestionExponentFunction\", and \"InteractionFunction\" \
(default Automatic, meaning use the current global MFGraphs` definitions of V, alpha, and g). \
When Tolerance > 0, iteration stops early when the infinity-norm change in flow variables between consecutive steps falls below the given tolerance.";

IsNonLinearSolution::usage =
"IsNonLinearSolution[Eqs] extracts AssoNonCritical and checks equations and inequalities. The right and left hand sides of the nonlinear equations are shown with the sup-norm of the difference."

H::usage =
"H[xi, p, m, edge] is the Hamiltonian function for the edges.
edge is a directed edge from the graph."

V::usage =
"V[x, edge] is the along-edge potential term in the Hamiltonian. Lower values make locations on an edge more attractive to agents. The default is 0.";

alpha::usage =
"alpha[edge] is the edge-dependent congestion exponent used in the Hamiltonian. The default is 1.";

g::usage =
"g[m, edge] is the edge-dependent interaction term as a function of density m. The default is -1/m^2.";

WithHamiltonianFunctions::usage =
"WithHamiltonianFunctions[vFun, alphaFun, gFun, expr] temporarily overrides the public MFGraphs Hamiltonian ingredients V, alpha, and g while evaluating expr. Any argument may be Automatic to keep the current definition.";

U::usage =
"U[x, edge, Eqs, sol] computes the value function at position x on the given edge."

PrecomputeM::usage =
"PrecomputeM[jMin, jMax, edge, nPoints] precomputes M[j,x,edge] on a grid and returns
an InterpolatingFunction. This avoids per-point FindRoot calls during NIntegrate.
nPoints controls the grid resolution (default 50).";

FastIntegratedMass::usage =
"FastIntegratedMass[interpM, j, edge] computes IntegratedMass using a precomputed interpolation of M.
interpM should be the result of PrecomputeM.";

IsFeasible::usage =
"IsFeasible[result] returns True if a solver result is feasible, checking legacy \
\"Status\" or standardized \"Feasibility\" keys.";

Options[NonLinearSolver] = {
    "MaxIterations" -> 15,
    "Tolerance" -> 0,
    "ReturnShape" -> "Legacy",
    "PotentialFunction" -> Automatic,
    "CongestionExponentFunction" -> Automatic,
    "InteractionFunction" -> Automatic
};

SetAttributes[WithHamiltonianFunctions, HoldRest];

WithHamiltonianFunctions[vFun_:Automatic, alphaFun_:Automatic, gFun_:Automatic, expr_] :=
    Block[{V, alpha, g},
        V     = If[vFun     =!= Automatic, vFun,     Function[{x, edge}, 0]];
        alpha = If[alphaFun =!= Automatic, alphaFun, 1&];
        g     = If[gFun     =!= Automatic, gFun,     Function[{m, edge}, -1/m^2]];
        expr
    ];

(* --- Hamiltonian and related functions --- *)

(* Default congestion exponent *)
alpha[edge_] := 1

(* Default interaction potential *)
g[m_, edge_] := -1/m^2

(* Default baseline: zero potential on every edge. *)
V[x_, edge_] := 0

H = Function[{xi, p, m, edge}, p^2/(2 m^alpha[edge]) + V[xi, edge] - g[m, edge]];

Begin["`Private`"];

(* --- NonLinearSolver: main iterative solver --- *)

NonLinearSolver[Eqs_, OptionsPattern[]] :=
    Module[ {potentialFunction, congestionExponentFunction, interactionFunction},
        potentialFunction = OptionValue["PotentialFunction"];
        congestionExponentFunction = OptionValue["CongestionExponentFunction"];
        interactionFunction = OptionValue["InteractionFunction"];
        WithHamiltonianFunctions[
            potentialFunction,
            congestionExponentFunction,
            interactionFunction,
            Module[ {AssoCritical, PreEqs = Eqs, AssoNonCritical, NonCriticalList, js,
             MaxIter = OptionValue["MaxIterations"], tol = OptionValue["Tolerance"],
             returnShape = OptionValue["ReturnShape"], status, flowKeys, flowVals,
             resultKind, message, solution, comparisonData},
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
                (* Feasibility check on the non-linear solution *)
                If[AssoNonCritical === Null,
                    status = "Infeasible",
                    status = CheckFlowFeasibility[AssoNonCritical]
                ];
                If[returnShape === "Standard",
                    resultKind = If[AssoNonCritical === Null, "Failure", "Success"];
                    message = If[AssoNonCritical === Null, "NoSolution", None];
                    solution = If[AssociationQ[AssoNonCritical], AssoNonCritical, Missing["NotAvailable"]];
                    comparisonData = BuildSolverComparisonData[PreEqs, solution];
                    Join[
                        PreEqs,
                        MakeSolverResult[
                            "NonLinear",
                            resultKind,
                            status,
                            message,
                            solution,
                            Join[comparisonData, <|
                                "AssoCritical" -> AssoCritical,
                                "AssoNonCritical" -> AssoNonCritical,
                                "Status" -> status
                            |>]
                        ]
                    ],
                    Join[PreEqs, <|"AssoCritical" -> AssoCritical, "AssoNonCritical" -> AssoNonCritical, "Status" -> status|>]
                ]
            ]
        ]
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
        MFGPrint[Style["Max error for non-linear solution: ", Bold, Blue], Max[Abs[Flatten[{Newlhs-Newrhs}]]]];
        approx
    ];

(* --- IsNonLinearSolution: validation --- *)

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
        Print[styleBlue@"Max error for non-linear solution: ", Max[Abs[Flatten[{N[Nlhs/.assoc]-(Nrhs/.assoc)}]]]];
        N@assoc
    ];

(* 
   Density recovery (M[j, x, edge]): Numerically solves the Hamilton-Jacobi equation 
   to find the agent distribution m corresponding to a given flow j and location x.
   
   This implements the "Current Method" or "Mass-Root Discovery" step described 
   in Section 3 of the Al Saleh et al. paper. This step is numerically delicate; 
   it uses multiple attempts with varying starting guesses and precisions 
   to ensure convergence to a positive real mass value. 
*)
SolveMassRoot[j_?NumericQ, x_?NumericQ, edge_] :=
    Module[{attempts, logAttempts, rule, value, wp, guess, accuracyGoal, precisionGoal,
      xVal, jVal, lowerBound, upperBound},
        attempts = {
            {MachinePrecision, 1., 8, 8},
            {MachinePrecision, 0.25, 8, 8},
            {MachinePrecision, 4., 8, 8},
            {30, 1., 12, 12},
            {30, 0.25, 12, 12},
            {30, 4., 12, 12},
            {50, 1., 16, 16}
        };
        logAttempts = {
            {MachinePrecision, -2., 8, 8},
            {MachinePrecision, 0., 8, 8},
            {MachinePrecision, 2., 8, 8},
            {30, -2., 12, 12},
            {30, 0., 12, 12},
            {30, 2., 12, 12},
            {50, 0., 16, 16}
        };
        lowerBound = 10^-10;
        upperBound = 10^6;
        Do[
            {wp, guess, accuracyGoal, precisionGoal} = attempt;
            xVal = SetPrecision[x, wp];
            jVal = SetPrecision[j, wp];
            rule = Quiet[
                FindRoot[
                    SetPrecision[MFGraphs`H[xVal, -jVal m^(MFGraphs`alpha[edge] - 1), m, edge], wp],
                    {m, SetPrecision[guess, wp], SetPrecision[lowerBound, wp], SetPrecision[upperBound, wp]},
                    WorkingPrecision -> wp,
                    AccuracyGoal -> accuracyGoal,
                    PrecisionGoal -> precisionGoal,
                    MaxIterations -> 100
                ],
                {
                    FindRoot::precw,
                    FindRoot::lstol,
                    FindRoot::reged,
                    FindRoot::cvmit,
                    FindRoot::jsing,
                    FindRoot::nlnum,
                    FindRoot::frsec
                }
            ];
            If[MatchQ[rule, {_Rule..}],
                value = m /. rule;
                If[
                    NumericQ[value] &&
                    TrueQ[Im[N[value]] == 0] &&
                    TrueQ[Re[N[value]] > 0] &&
                    Not[PossibleZeroQ[N[value] - lowerBound]] &&
                    Not[PossibleZeroQ[N[value] - upperBound]],
                    Return[N[Re[value]]]
                ]
            ],
            {attempt, attempts}
        ];
        Do[
            {wp, guess, accuracyGoal, precisionGoal} = attempt;
            xVal = SetPrecision[x, wp];
            jVal = SetPrecision[j, wp];
            rule = Quiet[
                FindRoot[
                    SetPrecision[MFGraphs`H[xVal, -jVal Exp[y]^(MFGraphs`alpha[edge] - 1), Exp[y], edge], wp],
                    {y, SetPrecision[guess, wp]},
                    WorkingPrecision -> wp,
                    AccuracyGoal -> accuracyGoal,
                    PrecisionGoal -> precisionGoal,
                    MaxIterations -> 100
                ],
                {
                    FindRoot::precw,
                    FindRoot::lstol,
                    FindRoot::reged,
                    FindRoot::cvmit,
                    FindRoot::jsing,
                    FindRoot::nlnum,
                    FindRoot::frsec
                }
            ];
            If[MatchQ[rule, {_Rule..}],
                value = Exp[y /. rule];
                If[NumericQ[value] && TrueQ[Im[N[value]] == 0] && TrueQ[Re[N[value]] > 0],
                    Return[N[Re[value]]]
                ]
            ],
            {attempt, logAttempts}
        ];
        0.
    ];

U[x_?NumericQ , edge_, Eqs_Association, sol_] :=
    Module[ {jay, uT},
        jay = Eqs["SignedFlows"][List@@edge]/.sol;
        uT = u@@Reverse@(List@@edge);
        uT = uT/.sol;
        If[ PossibleZeroQ[jay],
            uT,
            uT - jay NIntegrate[M[jay, y, edge]^(MFGraphs`alpha[edge] - 1), {y, 0, x}]
        ]
    ]

M[j_?NumericQ, x_?NumericQ, edge_] :=
    If[ PossibleZeroQ[j],
        0.,
        SolveMassRoot[j, x, edge]
    ];

IntegratedMass[j_?NumericQ, edge_] :=
    If[ PossibleZeroQ[j],
        0,
        j NIntegrate[M[j, x, edge]^(MFGraphs`alpha[edge] - 1), {x, 0, 1}] // Quiet
    ];

Cost[j_, edge_] :=
    IntegratedMass[j, edge];

(* --- Interpolation-based M for faster Cost evaluation --- *)

PrecomputeM[jMin_?NumericQ, jMax_?NumericQ, edge_, nPoints_Integer:50] :=
  Module[{jVals, xVals, mGrid, flatData},
    jVals = Subdivide[jMin, jMax, nPoints];
    xVals = Subdivide[0., 1., nPoints];
    EnsureParallelKernels[];
    If[!TrueQ[$MFGraphsParallelReady],
      DistributeDefinitions["MFGraphs`", "MFGraphs`Private`"];
      $MFGraphsParallelReady = True
    ];
    mGrid = ParallelTable[
      If[PossibleZeroQ[jv], 0.,
        SolveMassRoot[jv, xv, edge]
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

FastIntegratedMass[interpM_, j_?NumericQ, edge_] :=
  If[PossibleZeroQ[j], 0,
    j NIntegrate[interpM[j, x]^(MFGraphs`alpha[edge]-1), {x, 0, 1}] // Quiet
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

(* --- Feasibility check --- *)

IsFeasible[result_Association] :=
    Lookup[result, "Feasibility", Lookup[result, "Status", "Unknown"]] === "Feasible";

End[];
