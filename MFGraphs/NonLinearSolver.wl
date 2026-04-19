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
Options: \"MaxIterations\" (default 15), \"Tolerance\" (default 0), \
\"PotentialFunction\", \"CongestionExponentFunction\", and \"InteractionFunction\" \
(default Automatic, meaning use the current global MFGraphs` definitions of V, alpha, and g). \
It always returns a standardized solver-result association containing solver metadata, \
comparison fields, convergence data, and the solver-specific payload keys. \
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

NormalizeEdgeFunction::usage =
"NormalizeEdgeFunction[spec, defaultVal] normalizes a user-provided per-edge specification \
(scalar, Association, Function, or Automatic) into a guaranteed Function[edge, ...]. \
Unlisted edges in an Association default to defaultVal (default 1).";

WithHamiltonianFunctions::usage =
"WithHamiltonianFunctions[vFun, alphaFun, gFun, expr] temporarily overrides the public MFGraphs Hamiltonian ingredients V, alpha, and g while evaluating expr. Any argument may be Automatic to keep the current definition. \
alphaFun is normalized via NormalizeEdgeFunction so that scalar, Association, and Function specs all work.";

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

ClassifyKKT::usage =
"ClassifyKKT[result] classifies a NonLinearSolver result against three KKT \
conditions and returns a four-key Association:
  \"KKTClass\"      -> \"Feasible\" | \"Borderline\" | \"Infeasible\"
  \"KKTReason\"     -> string naming the binding condition, or None when Feasible
  \"KKTMetrics\"    -> <|\"PrimalMinFlow\", \"KirchhoffResidual\", \"ConvergenceResidual\"|>
  \"KKTViolations\" -> list of violated condition names (subset of {\"Primal\", \
\"Kirchhoff\", \"Convergence\"})

The three conditions and their two-tier thresholds:
  Primal      - min flow value must be >= PrimalFeasibleThreshold (default 0); \
within PrimalBorderlineThreshold (default -10^-6) is Borderline; below that, Infeasible.
  Kirchhoff   - flow conservation residual ||Kj - b||_inf must be <= \
KirchhoffFeasibleThreshold (default 10^-6); up to KirchhoffBorderlineThreshold \
(default 10^-3) is Borderline; above that, Infeasible.
  Convergence - final flow-change residual from iteration must be <= \
ConvergenceFeasibleThreshold (default 10^-8); up to ConvergenceBorderlineThreshold \
(default 10^-4) is Borderline; above that, Infeasible. Missing residual (e.g. exact \
fixed-point reached in one step) is treated as Feasible.

Classification rule: a result is Infeasible if any condition exceeds its Borderline \
threshold; Borderline if no condition is Infeasible but at least one exceeds its \
Feasible threshold or a metric is uncomputable; Feasible otherwise.

ClassifyKKT is read-only: it never modifies the input result, and all metrics are \
extracted from data already present in the result envelope (no expensive recomputation).";

Options[ClassifyKKT] = {
    (* Primal: min(j[e]) thresholds. Values above PrimalFeasibleThreshold are fully feasible. *)
    "PrimalFeasibleThreshold"      -> 0,
    "PrimalBorderlineThreshold"    -> -10^-6,
    (* Kirchhoff: ||Kj - b||_inf thresholds. Values at or below Feasible are fully feasible. *)
    "KirchhoffFeasibleThreshold"   -> 10^-6,
    "KirchhoffBorderlineThreshold" -> 10^-3,
    (* Convergence: final flow-change residual between last two iterates. *)
    "ConvergenceFeasibleThreshold"   -> 10^-8,
    "ConvergenceBorderlineThreshold" -> 10^-4
};

Options[NonLinearSolver] = {
    "MaxIterations" -> 15,
    "Tolerance" -> 0,
    "PotentialFunction" -> Automatic,
    "CongestionExponentFunction" -> Automatic,
    "InteractionFunction" -> Automatic
};

NormalizeEdgeFunction[spec_, defaultVal_: 1] := Which[
    spec === Automatic,    With[{d = defaultVal}, Function[edge, d]],
    NumericQ[spec],        With[{v = spec}, Function[edge, v]],
    AssociationQ[spec],    With[{a = spec, d = defaultVal}, Function[edge, Lookup[a, Key[edge], d]]],
    True,                  spec
];

SetAttributes[WithHamiltonianFunctions, HoldRest];

WithHamiltonianFunctions[vFun_:Automatic, alphaFun_:Automatic, gFun_:Automatic, expr_] :=
    Block[{V, alpha, g},
        V     = If[vFun     =!= Automatic, vFun,     Function[{x, edge}, 0]];
        alpha = NormalizeEdgeFunction[alphaFun, 1];
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
             status, resultKind, message, solution, comparisonData, convergenceData,
             flowDelta = Missing["NotAvailable"], iterations, stopReason, solveTime,
             seedMethod, timingResult, residualHistory = {}, demotedDueToFeasibilityQ},
                {solveTime, timingResult} = AbsoluteTiming[
                    If[ KeyExistsQ[PreEqs, "AssoNonCritical"],
                        AssoNonCritical = PreEqs["AssoNonCritical"];
                        seedMethod = "ProvidedNonLinearSeed"
                        ,
                        If[KeyExistsQ[PreEqs, "AssoCritical"],
                            AssoNonCritical = PreEqs["AssoCritical"];
                            seedMethod = "ProvidedCriticalSeed"
                            ,
                            PreEqs = CriticalCongestionSolver[PreEqs];
                            AssoNonCritical = Lookup[PreEqs, "AssoCritical", Null];
                            seedMethod =
                                If[AssociationQ[AssoNonCritical],
                                    "CriticalCongestionSeed",
                                    "CriticalSeedUnavailable"
                                ]
                        ]
                    ];
                js = Lookup[PreEqs, "js", {}];
                NonCriticalList =
                    If[!AssociationQ[AssoNonCritical],
                        {Null},
                        If[ tol > 0 && js =!= {},
                            (* Tolerance-based stopping: compare consecutive flow vectors *)
                            NestWhileList[
                                nonLinearStep[PreEqs],
                                AssoNonCritical,
                                (Norm[N[Values[KeyTake[#2, js]] - Values[KeyTake[#1, js]]], Infinity] > tol)&,
                                2, MaxIter
                            ],
                            (* Default: exact fixed-point check (stops when consecutive results are identical) *)
                            FixedPointList[nonLinearStep[PreEqs], AssoNonCritical, MaxIter]
                        ]
                    ];
                ];
                MFGPrint["Iterated ", Length[NonCriticalList]-1, " times out of ", MaxIter];
                AssoCritical = Lookup[
                    PreEqs,
                    "AssoCritical",
                    If[Length[NonCriticalList] >= 2, NonCriticalList[[2]], First[NonCriticalList]]
                ];
                AssoNonCritical = NonCriticalList // Last;
                iterations = Max[Length[NonCriticalList] - 1, 0];
                If[iterations > 0 && js =!= {},
                    flowDelta = Quiet @ Check[
                        Norm[
                            N[Values[KeyTake[NonCriticalList[[-1]], js]] - Values[KeyTake[NonCriticalList[[-2]], js]]],
                            Infinity
                        ],
                        Missing["NotAvailable"]
                    ]
                ];
                residualHistory =
                    If[Length[NonCriticalList] >= 2 && js =!= {} &&
                        AllTrue[NonCriticalList, AssociationQ],
                        Quiet @ Table[
                            Check[
                                Norm[
                                    N[Values[KeyTake[NonCriticalList[[k]], js]]
                                      - Values[KeyTake[NonCriticalList[[k-1]], js]]],
                                    Infinity
                                ],
                                Missing["NotAvailable"]
                            ],
                            {k, 2, Length[NonCriticalList]}
                        ],
                        {}
                    ];
                (* Feasibility check on the non-linear solution *)
                If[AssoNonCritical === Null,
                    status = "Infeasible",
                    status = CheckFlowFeasibility[AssoNonCritical]
                ];
                stopReason =
                    Which[
                        AssoNonCritical === Null, "NoSolution",
                        tol > 0 && NumericQ[flowDelta] && flowDelta <= tol, "ToleranceMet",
                        tol > 0, "MaxIterationsReached",
                        iterations > 0 && NonCriticalList[[-1]] === NonCriticalList[[-2]], "FixedPointReached",
                        iterations >= MaxIter, "MaxIterationsReached",
                        True, "IterationStopped"
                    ];
                resultKind =
                    Which[
                        AssoNonCritical === Null, "Failure",
                        MemberQ[{"ToleranceMet", "FixedPointReached"}, stopReason], "Success",
                        True, "NonConverged"
                    ];
                demotedDueToFeasibilityQ = resultKind === "Success" && status =!= "Feasible";
                If[demotedDueToFeasibilityQ, resultKind = "NonConverged"];
                message =
                    Which[
                        resultKind === "Success", None,
                        demotedDueToFeasibilityQ, "FlowFeasibilityFailed",
                        True, stopReason
                    ];
                solution = If[AssociationQ[AssoNonCritical], AssoNonCritical, Missing["NotAvailable"]];
                comparisonData = BuildSolverComparisonData[PreEqs, solution];
                convergenceData = <|
                    "StopReason" -> stopReason,
                    "FinalResidual" -> flowDelta,
                    "ResidualHistory" -> residualHistory,
                    "Iterations" -> iterations,
                    "SolveTime" -> solveTime,
                    "SeedMethod" -> seedMethod
                |>;
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
                            "Status" -> status,
                            "Convergence" -> convergenceData
                        |>]
                    ]
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

(* --- KKT gate --- *)

(* Match flow keys by head symbol name so both MFGraphs`j[...] and Global`j[...] are accepted. *)
kktFlowKeyQ[key_] := MatchQ[key, _Symbol[__]] && SymbolName[Head[key]] === "j";

(* Extract the three KKT metrics from an already-computed result envelope.
   All three read fields that NonLinearSolver writes unconditionally, so no
   expensive recomputation (NIntegrate, FindRoot, etc.) is triggered. *)
kktExtractMetrics[result_Association] :=
    Module[{solution, flowAssoc, flowVals, minFlow, kirchhoff, convergence},
        (* Primal: minimum value among all j[...] flow variables in the solution *)
        solution = Lookup[result, "AssoNonCritical",
                       Lookup[result, "Solution", Missing["NotAvailable"]]];
        minFlow =
            If[AssociationQ[solution],
                (
                    flowAssoc = KeySelect[solution, kktFlowKeyQ];
                    flowVals = Values[flowAssoc];
                    Which[
                        flowVals === {}, Missing["NotAvailable"],
                        !VectorQ[flowVals, NumericQ], Missing["NotAvailable"],
                        True, N @ Min[flowVals]
                    ]
                ),
                Missing["NotAvailable"]
            ];
        (* Conservation: Kirchhoff residual written by BuildSolverComparisonData *)
        kirchhoff = Lookup[result, "KirchhoffResidual", Missing["NotAvailable"]];
        (* Convergence: final flow-change between last two iterates, written by NonLinearSolver *)
        convergence =
            Lookup[
                Lookup[result, "Convergence", <||>],
                "FinalResidual",
                Missing["NotAvailable"]
            ];
        <|
            "PrimalMinFlow"       -> minFlow,
            "KirchhoffResidual"   -> kirchhoff,
            "ConvergenceResidual" -> convergence
        |>
    ];

(* Classify a single metric value against its two-tier thresholds.
   Returns "Feasible", "Borderline", or "Infeasible", with "Borderline" for Missing. *)
kktClassifyMetric[value_, feasibleTol_, borderlineTol_, lowerIsBetter_:True] :=
    Which[
        MissingQ[value] || !NumericQ[value], "Borderline",
        lowerIsBetter,
            Which[
                value <= feasibleTol,   "Feasible",
                value <= borderlineTol, "Borderline",
                True,                   "Infeasible"
            ],
        (* higher is better (primal: higher min-flow is better) *)
        True,
            Which[
                value >= feasibleTol,   "Feasible",
                value >= borderlineTol, "Borderline",
                True,                   "Infeasible"
            ]
    ];

ClassifyKKT[result_Association, opts:OptionsPattern[]] :=
    Module[{
        primalFeasTol, primalBorderTol,
        kirchhoffFeasTol, kirchhoffBorderTol,
        convFeasTol, convBorderTol,
        metrics,
        primalClass, kirchhoffClass, convClass,
        violations, kktClass, reason, anyMetricMissingQ
    },
        primalFeasTol      = OptionValue["PrimalFeasibleThreshold"];
        primalBorderTol    = OptionValue["PrimalBorderlineThreshold"];
        kirchhoffFeasTol   = OptionValue["KirchhoffFeasibleThreshold"];
        kirchhoffBorderTol = OptionValue["KirchhoffBorderlineThreshold"];
        convFeasTol        = OptionValue["ConvergenceFeasibleThreshold"];
        convBorderTol      = OptionValue["ConvergenceBorderlineThreshold"];

        metrics = kktExtractMetrics[result];
        anyMetricMissingQ = AnyTrue[Values[metrics], MissingQ];

        primalClass    = kktClassifyMetric[metrics["PrimalMinFlow"],
                             primalFeasTol, primalBorderTol, False (* higher is better *)];
        kirchhoffClass = kktClassifyMetric[metrics["KirchhoffResidual"],
                             kirchhoffFeasTol, kirchhoffBorderTol, True (* lower is better *)];
        (* Missing convergence residual means exact fixed-point in FixedPointList — treat as Feasible *)
        convClass =
            If[MissingQ[metrics["ConvergenceResidual"]],
                "Feasible",
                kktClassifyMetric[metrics["ConvergenceResidual"],
                    convFeasTol, convBorderTol, True (* lower is better *)]
            ];

        violations = Join[
            If[primalClass    =!= "Feasible", {"Primal"},      {}],
            If[kirchhoffClass =!= "Feasible", {"Kirchhoff"},   {}],
            If[convClass      =!= "Feasible", {"Convergence"}, {}]
        ];

        kktClass =
            Which[
                MemberQ[{primalClass, kirchhoffClass, convClass}, "Infeasible"], "Infeasible",
                violations =!= {},                                               "Borderline",
                True,                                                            "Feasible"
            ];

        reason =
            Which[
                kktClass === "Feasible",            None,
                primalClass    === "Infeasible",    "NegativeFlow",
                kirchhoffClass === "Infeasible",    "KirchhoffViolation",
                convClass      === "Infeasible",    "NotConverged",
                anyMetricMissingQ,                  "MetricsUnavailable",
                primalClass    === "Borderline",    "SmallNegativeFlow",
                kirchhoffClass === "Borderline",    "SmallKirchhoffViolation",
                convClass      === "Borderline",    "SlowConvergence",
                True,                               "MetricsUnavailable"
            ];

        <|
            "KKTClass"      -> kktClass,
            "KKTReason"     -> reason,
            "KKTMetrics"    -> metrics,
            "KKTViolations" -> violations
        |>
    ];

End[];
