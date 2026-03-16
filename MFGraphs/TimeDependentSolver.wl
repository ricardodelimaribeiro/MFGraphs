(* Wolfram Language package *)
(* Time-dependent MFG solver: backward-forward sweep on a time-discretized grid *)

(* --- Public API declarations --- *)

TimeDependentSolver::usage =
"TimeDependentSolver[data] solves the time-dependent MFG problem on a network.
The data Association must include \"Time Horizon\" and optionally
\"Time Steps\", \"Initial Mass Distribution\", \"Terminal Cost Function\",
\"Time Dependent Entrance Flows\", and \"Time Dependent Switching Costs\".
Options: \"MaxOuterIterations\" (default 20), \"Tolerance\" (default 10^-6),
\"SpatialSolverIterations\" (default 0; when > 0, runs nonlinear iterations per step),
\"ReturnShape\" (default \"Legacy\"; use \"Standard\" for normalized output).";

TimeDependentSolver::nottimedep =
"Data does not contain time-dependent keys. Use \"Time Horizon\" to specify the time domain.";
TimeDependentSolver::validation =
"Data validation failed: `1`.";

IsTimeDependentQ::usage =
"IsTimeDependentQ[data] returns True if the data Association contains time-dependent keys
(specifically, \"Time Horizon\").";

ValidateTimeDependentData::usage =
"ValidateTimeDependentData[data] checks that the required time-dependent keys are present
and consistent. Returns True or a failure message string.";

Options[TimeDependentSolver] = {
    "MaxOuterIterations" -> 20,
    "Tolerance" -> 10^-6,
    "SpatialSolverIterations" -> 0,
    "ReturnShape" -> "Legacy"
};

Begin["`Private`"];

(* ============================================================ *)
(* Detection and validation                                      *)
(* ============================================================ *)

IsTimeDependentQ[data_Association] :=
    KeyExistsQ[data, "Time Horizon"]

IsTimeDependentQ[_] := False

ValidateTimeDependentData[data_Association] :=
    Module[{T, Nt},
        If[!KeyExistsQ[data, "Time Horizon"],
            Return["Missing key: \"Time Horizon\"", Module]];
        T = data["Time Horizon"];
        If[!NumericQ[T] || T <= 0,
            Return["\"Time Horizon\" must be a positive number", Module]];
        Nt = Lookup[data, "Time Steps", 10];
        If[!IntegerQ[Nt] || Nt < 1,
            Return["\"Time Steps\" must be a positive integer", Module]];
        If[!KeyExistsQ[data, "Vertices List"],
            Return["Missing required network key: \"Vertices List\"", Module]];
        If[!KeyExistsQ[data, "Adjacency Matrix"],
            Return["Missing required network key: \"Adjacency Matrix\"", Module]];
        True
    ]

(* ============================================================ *)
(* Bridge function: rebuild d2e for a specific time step         *)
(* ============================================================ *)

(* Rebuilds the d2e Association with new entrance flows and exit costs,
   then drops cached preprocessing to force MFGPreprocessing to recompute. *)

rebuildD2EForTimeStep[d2e_Association, entranceFlowValues_List, exitCostValues_List] :=
    Module[{modified = d2e,
            inAuxEntryPairs, outAuxExitPairs, auxExitVertices,
            newEntryDataAssociation, newEqEntryIn, newRuleExitValues,
            newExitCosts},

        (* Structural info from d2e *)
        inAuxEntryPairs = List @@@ Lookup[d2e, "entryEdges", {}];
        outAuxExitPairs = List @@@ Lookup[d2e, "exitEdges", {}];
        auxExitVertices = Lookup[d2e, "auxExitVertices", {}];

        (* Rebuild entry flow data *)
        newEntryDataAssociation = RoundValues @ AssociationThread[
            inAuxEntryPairs, entranceFlowValues];
        newEqEntryIn = (j @@ # == newEntryDataAssociation[#]) & /@ inAuxEntryPairs;

        (* Rebuild exit cost data *)
        newRuleExitValues = AssociationThread[
            u @@@ (Reverse /@ outAuxExitPairs), exitCostValues];
        newRuleExitValues = Join[newRuleExitValues,
            AssociationThread[u @@@ outAuxExitPairs, exitCostValues]];
        newExitCosts = AssociationThread[auxExitVertices, exitCostValues];

        (* Update d2e keys *)
        modified["EntryDataAssociation"] = newEntryDataAssociation;
        modified["EqEntryIn"] = newEqEntryIn;
        modified["RuleExitValues"] = newRuleExitValues;
        modified["ExitCosts"] = newExitCosts;

        (* Drop cached preprocessing to force recomputation *)
        KeyDrop[modified, {"InitRules", "NewSystem"}]
    ]

(* ============================================================ *)
(* Helpers: extract boundary values from a solution              *)
(* ============================================================ *)

(* Extract exit vertex values from a solved time-step association *)
extractExitCostValues[d2e_Association, solution_Association] :=
    Module[{outAuxExitPairs, exitKeys},
        outAuxExitPairs = List @@@ Lookup[d2e, "exitEdges", {}];
        exitKeys = u @@@ outAuxExitPairs;
        Lookup[solution, exitKeys, 0]
    ]

(* Get entrance flow values at time t from the data *)
getEntranceFlowsAtTime[data_Association, d2e_Association, t_] :=
    Module[{entranceFlowsFn, entranceFlowsAtT, entryVerticesFlows},
        entranceFlowsFn = Lookup[data, "Time Dependent Entrance Flows", None];
        If[entranceFlowsFn === None,
            (* Fall back to static entrance flows *)
            entryVerticesFlows = Lookup[data, "Entrance Vertices and Flows", {}];
            Last /@ entryVerticesFlows,
            (* Evaluate the time-dependent function *)
            entranceFlowsAtT = entranceFlowsFn[t];
            Last /@ entranceFlowsAtT
        ]
    ]

(* Get exit cost values at terminal time from the data *)
getTerminalCostValues[data_Association, d2e_Association] :=
    Module[{terminalCostFn, exitVerticesCosts, exitVertices},
        terminalCostFn = Lookup[data, "Terminal Cost Function", None];
        exitVerticesCosts = Lookup[data, "Exit Vertices and Terminal Costs", {}];
        If[terminalCostFn === None,
            (* Fall back to static terminal costs *)
            Last /@ exitVerticesCosts,
            (* Evaluate the terminal cost function at exit vertices *)
            (* For vertex-level terminal costs, evaluate at x=0 on the exit edge *)
            exitVertices = First /@ exitVerticesCosts;
            Table[terminalCostFn[0, {exitVertices[[i]]}], {i, Length[exitVertices]}]
        ]
    ]

(* ============================================================ *)
(* Spatial solve at a single time step                           *)
(* ============================================================ *)

solveAtTimeStep[d2e_Association, entranceFlowValues_List,
    exitCostValues_List, approxJs_Association, spatialIter_Integer] :=
    Module[{d2eStep, preStep, solution, js, time, temp},

        (* Rebuild d2e with time-step boundary data *)
        d2eStep = rebuildD2EForTimeStep[d2e, entranceFlowValues, exitCostValues];

        (* Preprocess (builds InitRules and NewSystem) *)
        ClearSolveCache[];
        temp = MFGPrintTemporary["  Preprocessing time step..."];
        {time, preStep} = AbsoluteTiming @ MFGPreprocessing[d2eStep];
        NotebookDelete[temp];
        MFGPrint["  Preprocessing took ", time, " seconds."];

        (* Solve the spatial system *)
        solution = MFGSystemSolver[preStep][approxJs];

        (* Optional nonlinear iterations *)
        If[spatialIter > 0 && AssociationQ[solution],
            js = Lookup[preStep, "js", {}];
            Do[
                ClearSolveCache[];
                solution = MFGSystemSolver[preStep][KeyTake[solution, js]],
                {spatialIter}
            ]
        ];

        solution
    ]

(* ============================================================ *)
(* Backward pass: solve HJB from t_N to t_0                     *)
(* ============================================================ *)

backwardPass[d2e_Association, data_Association, timeGrid_List,
    massField_List, spatialIter_Integer] :=
    Module[{Nt = Length[timeGrid] - 1, valueField, flowField,
            exitCostValues, entranceFlowValues, solution, js, approxJs,
            k, time},

        js = Lookup[d2e, "js", {}];
        valueField = ConstantArray[<||>, Nt + 1];
        flowField = ConstantArray[<||>, Nt + 1];

        (* Terminal step: solve at t_N with original terminal costs *)
        MFGPrint["Backward pass: step ", Nt, "/", Nt, " (t = ", timeGrid[[Nt + 1]], ")"];
        exitCostValues = getTerminalCostValues[data, d2e];
        entranceFlowValues = getEntranceFlowsAtTime[data, d2e, timeGrid[[Nt + 1]]];
        approxJs = AssociationThread[js, 0 js];

        solution = solveAtTimeStep[d2e, entranceFlowValues, exitCostValues,
            approxJs, spatialIter];

        If[AssociationQ[solution],
            valueField[[Nt + 1]] = solution;
            flowField[[Nt + 1]] = KeyTake[solution, js],
            MFGPrint["Warning: terminal time step failed to solve."];
            valueField[[Nt + 1]] = <||>;
            flowField[[Nt + 1]] = AssociationThread[js, 0 js]
        ];

        (* Sweep backward: t_{N-1}, ..., t_0 *)
        Do[
            MFGPrint["Backward pass: step ", k, "/", Nt, " (t = ", timeGrid[[k + 1]], ")"];

            (* Exit costs come from the solution at the next time step *)
            exitCostValues = extractExitCostValues[d2e, valueField[[k + 2]]];
            entranceFlowValues = getEntranceFlowsAtTime[data, d2e, timeGrid[[k + 1]]];

            (* Use flows from previous iteration's mass field as approximation *)
            approxJs = If[AssociationQ[massField[[k + 1]]],
                massField[[k + 1]],
                AssociationThread[js, 0 js]
            ];

            solution = solveAtTimeStep[d2e, entranceFlowValues, exitCostValues,
                approxJs, spatialIter];

            If[AssociationQ[solution],
                valueField[[k + 1]] = solution;
                flowField[[k + 1]] = KeyTake[solution, js],
                MFGPrint["Warning: time step ", k, " failed to solve."];
                valueField[[k + 1]] = <||>;
                flowField[[k + 1]] = AssociationThread[js, 0 js]
            ],

            {k, Nt - 1, 0, -1}
        ];

        {valueField, flowField}
    ]

(* ============================================================ *)
(* Forward pass: compute mass from flows                         *)
(* ============================================================ *)

forwardPass[d2e_Association, flowField_List, timeGrid_List] :=
    Module[{Nt = Length[timeGrid] - 1, massField, js},
        js = Lookup[d2e, "js", {}];

        (* For Phase 1: mass field stores the flow approximations.
           In a future phase, this will compute actual mass distributions
           via M[j, x, edge] and propagate via the FPK equation. *)
        massField = Table[
            If[AssociationQ[flowField[[k]]],
                flowField[[k]],
                AssociationThread[js, 0 js]
            ],
            {k, 1, Nt + 1}
        ];

        massField
    ]

(* ============================================================ *)
(* Result assembly                                               *)
(* ============================================================ *)

assembleResult[timeGrid_List, valueField_List, flowField_List,
    massField_List, iterations_Integer, residual_, converged_,
    returnShape_String] :=
    Module[{solution, status, resultKind, message,
            timeAssocValue, timeAssocFlow, timeAssocMass},

        (* Build time-indexed associations *)
        timeAssocValue = AssociationThread[timeGrid, valueField];
        timeAssocFlow = AssociationThread[timeGrid, flowField];
        timeAssocMass = AssociationThread[timeGrid, massField];

        solution = <|
            "TimeGrid" -> timeGrid,
            "ValueFunction" -> timeAssocValue,
            "FlowField" -> timeAssocFlow,
            "MassDistribution" -> timeAssocMass,
            "Convergence" -> <|
                "Iterations" -> iterations,
                "Residual" -> residual,
                "Converged" -> converged
            |>
        |>;

        status = If[converged, "Feasible", "Infeasible"];
        resultKind = If[converged, "Success", "Failure"];
        message = If[converged, None,
            "Did not converge in " <> ToString[iterations] <> " iterations"];

        If[returnShape === "Standard",
            MakeSolverResult[
                "TimeDependentMFG",
                resultKind,
                status,
                message,
                solution,
                <|"AssoTimeDependentMFG" -> solution, "Status" -> status|>
            ],
            <|"AssoTimeDependentMFG" -> solution, "Status" -> status|>
        ]
    ]

(* ============================================================ *)
(* Main solver: outer backward-forward iteration                 *)
(* ============================================================ *)

TimeDependentSolver[data_Association, opts:OptionsPattern[]] :=
    Module[{d2e, T, Nt, dt, timeGrid,
            maxIter, tol, spatialIter, returnShape,
            massField, valueField, flowField, prevFlowField,
            residual = Infinity, converged = False, iteration = 0,
            js, validation, time},

        (* Validate *)
        If[!IsTimeDependentQ[data],
            Message[TimeDependentSolver::nottimedep];
            Return[$Failed, Module]
        ];
        validation = ValidateTimeDependentData[data];
        If[validation =!= True,
            Message[TimeDependentSolver::validation, validation];
            Return[$Failed, Module]
        ];

        (* Extract options *)
        maxIter = OptionValue["MaxOuterIterations"];
        tol = OptionValue["Tolerance"];
        spatialIter = OptionValue["SpatialSolverIterations"];
        returnShape = OptionValue["ReturnShape"];

        (* Time discretization *)
        T = data["Time Horizon"];
        Nt = Lookup[data, "Time Steps", 10];
        dt = N[T / Nt];
        timeGrid = N @ Subdivide[0, T, Nt];

        (* Build spatial structure once *)
        MFGPrint["Building spatial structure..."];
        {time, d2e} = AbsoluteTiming @ DataToEquations[data];
        MFGPrint["DataToEquations took ", time, " seconds."];

        js = Lookup[d2e, "js", {}];

        (* Initialize: zero flows at all time steps *)
        massField = Table[AssociationThread[js, 0 js], {Nt + 1}];

        (* Outer backward-forward loop *)
        Do[
            MFGPrint["\n=== Outer iteration ", iteration, " / ", maxIter, " ==="];

            (* Backward pass: solve for value function *)
            {time, {valueField, flowField}} = AbsoluteTiming @
                backwardPass[d2e, data, timeGrid, massField, spatialIter];
            MFGPrint["Backward pass took ", time, " seconds."];

            (* Convergence check: compare flow fields *)
            If[iteration > 1,
                residual = Max @ Table[
                    If[AssociationQ[flowField[[k]]] && AssociationQ[prevFlowField[[k]]],
                        Module[{diff = Values[flowField[[k]]] - Values[prevFlowField[[k]]]},
                            If[AllTrue[diff, NumericQ],
                                Max[Abs[diff]],
                                Infinity
                            ]
                        ],
                        Infinity
                    ],
                    {k, 1, Nt + 1}
                ];
                MFGPrint["Flow residual: ", residual];

                If[residual < tol,
                    converged = True;
                    MFGPrint["Converged!"];
                    Break[]
                ]
            ];

            prevFlowField = flowField;

            (* Forward pass: update mass field from flows *)
            massField = forwardPass[d2e, flowField, timeGrid];
            MFGPrint["Forward pass complete."],

            {iteration, 1, maxIter}
        ];

        (* Assemble and return result *)
        assembleResult[timeGrid, valueField, flowField, massField,
            iteration, residual, converged, returnShape]
    ]

End[];
