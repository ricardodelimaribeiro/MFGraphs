(* Wolfram Language package *)
(*
   Solvers: Critical-congestion solver suite for MFGraphs.

   Extracted from DataToEquations.wl as part of the Phase 3 architectural
   refactor. Contains all solver-layer code: DirectCriticalSolver,
   CriticalCongestionSolver, numeric backends, IsCriticalSolution, and
   all supporting helpers.

   Load order dependency: DataToEquations.wl and DNFReduce.wl must be
   loaded before this file (enforced by MFGraphs.wl).
*)

BeginPackage["MFGraphs`"];

Begin["`Private`"];

(* --- Solver tolerance constant --- *)
(* Single source of truth for numeric tolerance across all critical-congestion
   solver paths: JFirst backend, coupled numeric backend, DirectCriticalSolver,
   and the IsCriticalSolution validation gate. The CriticalCongestionSolver
   "ValidationTolerance" option overrides this when set explicitly. *)
$CriticalSolverTolerance = 10^-6;

GraphDistanceHeuristicSafeQ[Eqs_Association, tol_:10^-10] :=
    Module[{exitCosts, exitCostSpread, switchingValues, nonzeroSwitchingQ},
        exitCosts = Last /@ Lookup[Eqs, "Exit Vertices and Terminal Costs", {}];
        exitCostSpread =
            Which[
                Length[exitCosts] <= 1, 0.,
                VectorQ[exitCosts, NumericQ], N @ Max[Abs[N @ exitCosts - First[N @ exitCosts]]],
                True, Infinity
            ];
        switchingValues = Values @ Lookup[Eqs, "SwitchingCosts", <||>];
        nonzeroSwitchingQ =
            AnyTrue[
                switchingValues,
                Function[val,
                    Which[
                        NumericQ[val], !PossibleZeroQ[val],
                        True, True
                    ]
                ]
            ];
        exitCostSpread <= tol && !nonzeroSwitchingQ
    ];

(* --- Graph-distance-based Or-resolution --- *)
(* Resolves flow complementarity conditions (j[a,b]==0 || j[b,a]==0)
   using shortest-path distances from exit vertices. The vertex farther
   from the exit sends flow toward the closer vertex, setting the
   backward flow to zero. This eliminates most Or-conditions before
   DNFReduce, reducing exponential branching to O(V+E) graph traversal. *)

ResolveOneOr[Or[lhs_, rhs_], distToExit_Association] :=
    Which[
        (* Edge flow: j[a,b]==0 || j[b,a]==0 *)
        MatchQ[{lhs, rhs}, {_j == 0, _j == 0}] &&
            Length[lhs[[1]]] == 2 && Length[rhs[[1]]] == 2 &&
            lhs[[1, 1]] === rhs[[1, 2]] && lhs[[1, 2]] === rhs[[1, 1]],
        With[{da = Lookup[distToExit, lhs[[1, 1]], Infinity],
              db = Lookup[distToExit, lhs[[1, 2]], Infinity]},
            Which[da > db, rhs, da < db, lhs, True, $Failed]
        ],
        (* Transition flow: j[a,b,c]==0 || j[c,b,a]==0 *)
        MatchQ[{lhs, rhs}, {_j == 0, _j == 0}] &&
            Length[lhs[[1]]] == 3 && Length[rhs[[1]]] == 3 &&
            lhs[[1, 1]] === rhs[[1, 3]] && lhs[[1, 3]] === rhs[[1, 1]] &&
            lhs[[1, 2]] === rhs[[1, 2]],
        With[{da = Lookup[distToExit, lhs[[1, 1]], Infinity],
              dc = Lookup[distToExit, lhs[[1, 3]], Infinity]},
            Which[da > dc, rhs, da < dc, lhs, True, $Failed]
        ],
        True, $Failed
    ]

ResolveOneOr[_, _] := $Failed

ResolveOrByGraphDistance[Eqs_Association, triple:{_, _, _}] :=
    Module[{graph, exitVerts, distToExit, orTerms, resolved = {},
            unresolved = {}, newEE, newOR},
        If[triple[[3]] === True, Return[triple, Module]];
        If[!GraphDistanceHeuristicSafeQ[Eqs],
            MFGPrint["ResolveOrByGraphDistance: skipped (heterogeneous exit/switching costs)."];
            Return[triple, Module]
        ];
        graph = Lookup[Eqs, "auxiliaryGraph", None];
        exitVerts = Lookup[Eqs, "auxExitVertices", {}];
        If[graph === None || exitVerts === {},
            MFGPrint["ResolveOrByGraphDistance: no graph data available"];
            Return[triple, Module]
        ];
        distToExit = Association @ Table[
            v -> Min[GraphDistance[graph, v, #]& /@ exitVerts],
            {v, VertexList[graph]}
        ];
        orTerms = If[Head[triple[[3]]] === And,
            List @@ triple[[3]], {triple[[3]]}];
        Do[
            With[{res = ResolveOneOr[term, distToExit]},
                If[res =!= $Failed,
                    AppendTo[resolved, res],
                    AppendTo[unresolved, term]
                ]
            ],
            {term, orTerms}
        ];
        MFGPrint["ResolveOrByGraphDistance: resolved ", Length[resolved],
            " of ", Length[orTerms], " Or-conditions"];
        newEE = If[resolved === {},
            triple[[1]],
            And @@ Join[
                If[triple[[1]] === True, {},
                    If[Head[triple[[1]]] === And, List @@ triple[[1]],
                        {triple[[1]]}]],
                resolved
            ]
        ];
        newOR = If[unresolved === {}, True, And @@ unresolved];
        {newEE, triple[[2]], newOR}
    ]

DirectCriticalSolver[Eqs_Association] :=
    Module[{graph, exitVerts, entryVerts, auxPairs, auxTriples,
            us, js, jts, distToExit,
            RuleEntryOut, RuleExitFlowsIn, RuleEntryValues, RuleExitValues,
            RuleBalanceGatheringFlows, EqGeneral, EqEntryIn,
            EqBalanceSplittingFlows, EqValueAuxiliaryEdges,
            AltFlows, AltTransitionFlows,
            zeroRules, eqSystem, allVars, sol, result},
        graph = Lookup[Eqs, "auxiliaryGraph", None];
        exitVerts = Lookup[Eqs, "auxExitVertices", {}];
        entryVerts = Lookup[Eqs, "auxEntryVertices", {}];
        us = Lookup[Eqs, "us", {}];
        js = Lookup[Eqs, "js", {}];
        jts = Lookup[Eqs, "jts", {}];
        auxPairs = Cases[us, u[a_, b_] :> {a, b}];
        auxTriples = Cases[jts, j[a_, b_, c_] :> {a, b, c}];
        RuleEntryOut = Lookup[Eqs, "RuleEntryOut", {}];
        RuleExitFlowsIn = Lookup[Eqs, "RuleExitFlowsIn", {}];
        RuleEntryValues = Lookup[Eqs, "RuleEntryValues", <||>];
        RuleExitValues = Lookup[Eqs, "RuleExitValues", <||>];
        RuleBalanceGatheringFlows = Lookup[Eqs, "RuleBalanceGatheringFlows", <||>];
        EqGeneral = Lookup[Eqs, "EqGeneral", True];
        EqEntryIn = Lookup[Eqs, "EqEntryIn", True];
        EqBalanceSplittingFlows = Lookup[Eqs, "EqBalanceSplittingFlows", True];
        EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", True];

        (* Phase 1: Compute vertex distances to exit *)
        distToExit = Association @ Table[
            v -> Min[GraphDistance[graph, v, #]& /@ exitVerts],
            {v, VertexList[graph]}
        ];

        (* Phase 2: Determine which flows are zero from graph distance.
           For edge {a,b}: if dist(a) < dist(b), then j[a,b]=0 (backward).
           For equidistant edges, both directions may carry flow.
           Transition flows are left free — determined by Kirchhoff. *)
        zeroRules = Association[];
        Do[With[{a = pair[[1]], b = pair[[2]]},
            With[{da = Lookup[distToExit, a, Infinity],
                  db = Lookup[distToExit, b, Infinity]},
                If[da < db, zeroRules[j[a, b]] = 0]
            ]], {pair, auxPairs}];
        MFGPrint["DirectCriticalSolver: ", Length[zeroRules],
            " flows set to zero by graph distance"];

        (* Phase 3: Build the full equation system.
           Substitute: Ncpc=0 (critical congestion at zero flow),
           known boundary rules, zero switching cost equalities
           (u[a,b]=u[c,b] at each vertex), and direction zeros. *)
        eqSystem = And @@ Flatten[{
            (* EqGeneral at Ncpc=0: j[a,b]-j[b,a]+u[a,b]-u[b,a]=0 *)
            List @@ (EqGeneral /. Thread[
                Cases[Keys[Lookup[Eqs, "costpluscurrents", <||>]], _] -> 0]),
            (* Entry flow *)
            If[Head[EqEntryIn] === And, List @@ EqEntryIn, {EqEntryIn}],
            (* Balance/splitting *)
            If[Head[EqBalanceSplittingFlows] === And,
                List @@ EqBalanceSplittingFlows,
                If[EqBalanceSplittingFlows === True, {}, {EqBalanceSplittingFlows}]],
            (* Value auxiliary edges *)
            If[Head[EqValueAuxiliaryEdges] === And,
                List @@ EqValueAuxiliaryEdges,
                If[EqValueAuxiliaryEdges === True, {}, {EqValueAuxiliaryEdges}]],
            (* Boundary values *)
            KeyValueMap[Equal, RuleExitValues],
            KeyValueMap[Equal, RuleEntryValues],
            (* Entry/exit flow rules *)
            (Equal @@@ Normal[RuleEntryOut]),
            (Equal @@@ Normal[RuleExitFlowsIn]),
            (* Gathering flow rules *)
            (Equal @@@ Normal[RuleBalanceGatheringFlows]),
            (* Direction zeros *)
            (Equal[#, 0]& /@ Keys[zeroRules]),
            (* Zero switching costs: u[a,b] = u[c,b] at each vertex *)
            (* Group by second element and equate all *)
            Module[{byVertex = GroupBy[auxPairs, Last]},
                Flatten @ KeyValueMap[
                    Function[{v, pairs},
                        If[Length[pairs] > 1,
                            (u @@ # == u @@ pairs[[1]])& /@ Rest[pairs],
                            {}
                        ]
                    ],
                    byVertex
                ]
            ]
        }];

        (* Phase 4: Solve with non-negativity on all flow variables *)
        allVars = Join[us, js, jts];
        With[{flowVars = Join[js, jts]},
            With[{ineqs = And @@ ((# >= 0 &) /@ flowVars)},
                With[{inst = Quiet @ FindInstance[eqSystem && ineqs, allVars, Reals]},
                    If[inst === {} || !MatchQ[inst, {{__Rule}, ___}],
                        MFGPrint["DirectCriticalSolver: FindInstance failed"];
                        Return[$Failed, Module]
                    ];
                    result = Association @ First @ inst;
                ]
            ]
        ];
        Do[If[!KeyExistsQ[result, var] && MatchQ[var, j[__]], result[var] = 0], {var, allVars}];
        (* Only resolve transitive chains if there are symbolic values *)
        If[!AllTrue[Values[result], NumericQ],
            result = Expand /@ FixedPoint[Function[r, ReplaceAll[r] /@ r], result, 10]
        ];
        result = Join[KeyTake[result, us], KeyTake[result, js], KeyTake[result, jts]];
        result
    ]

(* --- u-variable equivalence reduction (zero-switching-cost preprocessing) --- *)

(* BuildUEquivalenceReduction: when switching costs are zero at a vertex,
   all u[a,v] for that vertex v are equal.  We treat these equalities as
   undirected edges in a graph over u-variables and collapse each connected
   component to a single representative.  If a boundary value is known for
   any member of the class, we prefer that member as representative so the
   substitution injects the numeric constant directly. *)
BuildUEquivalenceReduction[Eqs_Association] :=
    Module[{us, auxPairs, switchingCosts, byVertex, edges, graph, components,
        boundaryVars, representative, rules, liftRules, nEliminated},
        us = Lookup[Eqs, "us", {}];
        auxPairs = Cases[us, u[a_, b_] :> {a, b}];
        switchingCosts = Lookup[Eqs, "SwitchingCosts", <||>];

        (* Group auxiliary pairs by destination vertex *)
        byVertex = GroupBy[auxPairs, Last];

        (* Build undirected equality edges: u[a,v] == u[c,v] when S(a,v,c)=0 *)
        edges = Flatten @ KeyValueMap[
            Function[{v, pairs},
                If[Length[pairs] <= 1, {},
                    Select[
                        Flatten[Table[
                            UndirectedEdge[u @@ pairs[[i]], u @@ pairs[[j]]],
                            {i, Length[pairs]}, {j, i + 1, Length[pairs]}
                        ]],
                        (* Both directions of switching must be zero *)
                        With[{p1 = List @@ #[[1]], p2 = List @@ #[[2]]},
                            Lookup[switchingCosts, {p1[[1]], v, p2[[1]]}, Infinity] === 0 &&
                            Lookup[switchingCosts, {p2[[1]], v, p1[[1]]}, Infinity] === 0
                        ] &
                    ]
                ]
            ],
            byVertex
        ];

        If[edges === {},
            Return[<|
                "SubstitutionRules" -> {},
                "LiftRules" -> {},
                "ReducedUVars" -> us,
                "EliminatedCount" -> 0,
                "ClassCount" -> Length[us],
                "Classes" -> (List /@ us)
            |>, Module]
        ];

        graph = Graph[us, edges];
        components = ConnectedComponents[graph];

        (* Identify boundary-valued u-variables for representative preference *)
        boundaryVars = Join[
            Keys[Lookup[Eqs, "RuleExitValues", <||>]],
            Keys[Lookup[Eqs, "RuleEntryValues", <||>]]
        ];

        (* For each class, pick representative: prefer boundary vars, else first *)
        rules = {};
        liftRules = {};
        Do[
            representative = With[{bnd = Select[class, MemberQ[boundaryVars, #] &]},
                If[bnd =!= {}, First[bnd], First[class]]
            ];
            Do[
                If[var =!= representative,
                    AppendTo[rules, var -> representative];
                    AppendTo[liftRules, var -> representative];
                ],
                {var, class}
            ],
            {class, components}
        ];

        nEliminated = Length[rules];
        <|
            "SubstitutionRules" -> rules,
            "LiftRules" -> liftRules,
            "ReducedUVars" -> Complement[us, Keys[Association[rules]]],
            "EliminatedCount" -> nEliminated,
            "ClassCount" -> Length[components],
            "Classes" -> components
        |>
    ];

(* --- Critical numeric backend (internal, guarded) --- *)

$CriticalNumericBackendEnabled = True;

CriticalNumericBackendRequestedQ[Eqs_Association] :=
    Module[{mode},
        mode = Lookup[Eqs, "CriticalNumericBackendMode", Automatic];
        Which[
            mode === True || mode === "Force", True,
            mode === False || mode === "Disabled", False,
            mode === Automatic || mode === "Auto", TrueQ[$CriticalNumericBackendEnabled],
            True, False
        ]
    ];

CriticalNumericBackendEligibleQ[Eqs_Association] :=
    Module[{numericState, switchingVals, entryVals, exitVals, km, b, vars},
        numericState = Lookup[Eqs, "NumericState", Missing["NotAvailable"]];
        switchingVals = Values[Lookup[Eqs, "SwitchingCosts", <||>]];
        entryVals = Last /@ Lookup[Eqs, "Entrance Vertices and Flows", {}];
        exitVals = Last /@ Lookup[Eqs, "Exit Vertices and Terminal Costs", {}];
        km = Lookup[numericState, "KirchhoffKM", Missing["NotAvailable"]];
        b = Lookup[numericState, "KirchhoffB", Missing["NotAvailable"]];
        vars = Lookup[numericState, "KirchhoffVariables", {}];
        AssociationQ[numericState] &&
        AllTrue[switchingVals, NumericQ[#] && PossibleZeroQ[#] &] &&
        Lookup[Eqs, "auxiliaryGraph", None] =!= None &&
        AllTrue[entryVals, NumericQ] &&
        AllTrue[exitVals, NumericQ] &&
        Head[km] === SparseArray &&
        VectorQ[b, NumericQ] &&
        ListQ[vars] &&
        vars =!= {}
    ];

BuildCriticalLinearConstraints[Eqs_Association] :=
    Module[{graph, exitVerts, us, js, jts, auxPairs, distToExit, zeroRules, pair,
        RuleEntryOut, RuleExitFlowsIn, RuleEntryValues, RuleExitValues,
        RuleBalanceGatheringFlows, EqGeneral, EqEntryIn,
        EqBalanceSplittingFlows, EqValueAuxiliaryEdges,
        eqList, flowVars, allVars, ineqList},
        graph = Lookup[Eqs, "auxiliaryGraph", None];
        exitVerts = Lookup[Eqs, "auxExitVertices", {}];
        us = Lookup[Eqs, "us", {}];
        js = Lookup[Eqs, "js", {}];
        jts = Lookup[Eqs, "jts", {}];
        auxPairs = Cases[us, u[a_, b_] :> {a, b}];
        RuleEntryOut = Lookup[Eqs, "RuleEntryOut", {}];
        RuleExitFlowsIn = Lookup[Eqs, "RuleExitFlowsIn", {}];
        RuleEntryValues = Lookup[Eqs, "RuleEntryValues", <||>];
        RuleExitValues = Lookup[Eqs, "RuleExitValues", <||>];
        RuleBalanceGatheringFlows = Lookup[Eqs, "RuleBalanceGatheringFlows", <||>];
        EqGeneral = Lookup[Eqs, "EqGeneral", True];
        EqEntryIn = Lookup[Eqs, "EqEntryIn", True];
        EqBalanceSplittingFlows = Lookup[Eqs, "EqBalanceSplittingFlows", True];
        EqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", True];

        distToExit = Lookup[
            Lookup[Eqs, "NumericState", <||>],
            "DistanceToExit",
            <||>
        ];
        If[!AssociationQ[distToExit] || distToExit === <||>,
            If[Head[graph] =!= Graph || exitVerts === {},
                Return[$Failed, Module]
            ];
            distToExit = Association @ Table[
                v -> Min[GraphDistance[graph, v, #] & /@ exitVerts],
                {v, VertexList[graph]}
            ];
        ];

        zeroRules = Association[];
        Do[
            pair = auxPairs[[k]];
            If[
                Lookup[distToExit, pair[[1]], Infinity] <
                    Lookup[distToExit, pair[[2]], Infinity],
                zeroRules[j @@ pair] = 0
            ],
            {k, Length[auxPairs]}
        ];

        eqList = Flatten[{
            List @@ (EqGeneral /. Thread[
                Keys[Lookup[Eqs, "costpluscurrents", <||>]] -> 0
            ]),
            If[Head[EqEntryIn] === And, List @@ EqEntryIn, If[EqEntryIn === True, {}, {EqEntryIn}]],
            If[Head[EqBalanceSplittingFlows] === And,
                List @@ EqBalanceSplittingFlows,
                If[EqBalanceSplittingFlows === True, {}, {EqBalanceSplittingFlows}]
            ],
            If[Head[EqValueAuxiliaryEdges] === And,
                List @@ EqValueAuxiliaryEdges,
                If[EqValueAuxiliaryEdges === True, {}, {EqValueAuxiliaryEdges}]
            ],
            KeyValueMap[Equal, RuleExitValues],
            KeyValueMap[Equal, RuleEntryValues],
            Equal @@@ Normal[RuleEntryOut],
            Equal @@@ Normal[RuleExitFlowsIn],
            Equal @@@ Normal[RuleBalanceGatheringFlows],
            (Equal[#, 0] & /@ Keys[zeroRules]),
            Module[{byVertex = GroupBy[auxPairs, Last]},
                Flatten @ KeyValueMap[
                    Function[{v, pairs},
                        If[
                            Length[pairs] > 1,
                            (u @@ # == u @@ pairs[[1]]) & /@ Rest[pairs],
                            {}
                        ]
                    ],
                    byVertex
                ]
            ]
        }];
        eqList = Cases[DeleteCases[eqList, True], _Equal];

        flowVars = Join[js, jts];
        allVars = Join[us, js, jts];
        ineqList = (# >= 0) & /@ flowVars;
        <|
            "Equalities" -> eqList,
            "Inequalities" -> ineqList,
            "FlowVariables" -> flowVars,
            "AllVariables" -> allVars
        |>
    ];

BuildLinearSystemFromEqualities[equalities_List, vars_List] :=
    Module[{exprs, placeholders, subExprs, coeffData, bVec, aMat},
        If[equalities === {} || vars === {},
            Return[$Failed, Module]
        ];
        exprs = equalities /. Equal[l_, r_] :> (l - r);
        placeholders = Table[Unique["lv"], {Length[vars]}];
        subExprs = exprs /. Thread[vars -> placeholders];
        coeffData = Quiet @ Check[
            CoefficientArrays[subExprs, placeholders],
            $Failed
        ];
        If[coeffData === $Failed || Length[coeffData] < 2,
            Return[$Failed, Module]
        ];
        bVec = Quiet @ Check[
            Developer`ToPackedArray @ N @ (-Normal[coeffData[[1]]]),
            $Failed
        ];
        aMat = Quiet @ Check[
            SparseArray[coeffData[[2]]],
            $Failed
        ];
        If[
            bVec === $Failed || aMat === $Failed || !VectorQ[bVec, NumericQ],
            $Failed,
            <|"A" -> aMat, "b" -> bVec|>
        ]
    ];

LinearSolveCandidate[aMat_SparseArray, bVec_List] :=
    Module[{dims, ls, x, ata, atb},
        dims = Dimensions[aMat];
        If[Length[dims] =!= 2 || dims[[2]] == 0,
            Return[$Failed, Module]
        ];
        If[dims[[1]] === dims[[2]],
            ls = Quiet @ Check[LinearSolve[N @ aMat], $Failed];
            If[ls === $Failed,
                Return[$Failed, Module]
            ];
            x = Quiet @ Check[N @ ls[N @ bVec], $Failed];
            If[VectorQ[x, NumericQ], Developer`ToPackedArray @ x, $Failed]
            ,
            ata = Quiet @ Check[N @ (Transpose[aMat].aMat), $Failed];
            atb = Quiet @ Check[N @ (Transpose[aMat].bVec), $Failed];
            If[ata === $Failed || atb === $Failed,
                Return[$Failed, Module]
            ];
            ls = Quiet @ Check[LinearSolve[ata], $Failed];
            If[ls === $Failed,
                x = Quiet @ Check[N @ LeastSquares[N @ Normal[aMat], N @ bVec], $Failed];
                Return[If[VectorQ[x, NumericQ], Developer`ToPackedArray @ x, $Failed], Module]
            ];
            x = Quiet @ Check[N @ ls[atb], $Failed];
            If[!VectorQ[x, NumericQ],
                x = Quiet @ Check[N @ LeastSquares[N @ Normal[aMat], N @ bVec], $Failed]
            ];
            If[VectorQ[x, NumericQ], Developer`ToPackedArray @ x, $Failed]
        ]
    ];

CriticalJFirstFailure[reason_String, details_Association : <||>] :=
    Failure["CriticalJFirstBackend", Join[<|"Reason" -> reason|>, details]];

CriticalJFirstFailureReason[result_] :=
    Which[
        MatchQ[result, Failure[_, _Association]],
            Lookup[result[[2]], "Reason", "JFirstFailed"],
        True,
            "JFirstFailed"
    ];

BuildCriticalJFirstObjectiveExpression[Eqs_Association, massVars_List] :=
    Module[{model, objective, objectiveVars},
        model = Quiet @ Check[BuildCriticalQuadraticObjective[Eqs], $Failed];
        If[!AssociationQ[model],
            Return[CriticalJFirstFailure["ObjectiveUnavailable"], Module]
        ];
        objective = Lookup[model, "Objective", Missing["NotAvailable"]];
        If[objective === Missing["NotAvailable"],
            Return[CriticalJFirstFailure["ObjectiveUnavailable"], Module]
        ];
        objectiveVars = Select[Variables[objective], MatchQ[#, _j] &];
        If[!SubsetQ[massVars, objectiveVars],
            Return[CriticalJFirstFailure["ObjectiveScopeMismatch"], Module]
        ];
        objective
    ];

CriticalConjunctionTerms[expr_] :=
    Which[
        TrueQ[expr], {},
        Head[expr] === And, List @@ expr,
        True, {expr}
    ];

RecoverCriticalFlowAssociation[Eqs_Association, flowBase_Association, tol_?NumericQ] :=
    Module[{entryInRules, flowRules, replacementRules, allFlowVars, resolve, resolved},
        entryInRules = Association @ Flatten[ToRules /@ Lookup[Eqs, "EqEntryIn", {}]];
        flowRules = Join[
            entryInRules,
            Lookup[Eqs, "RuleEntryOut", <||>],
            Lookup[Eqs, "RuleExitFlowsIn", <||>],
            Lookup[Eqs, "RuleBalanceGatheringFlows", <||>]
        ];
        replacementRules = Join[Normal[flowBase], Normal[flowRules]];
        resolve[expr_] := FixedPoint[
            Function[val, Quiet @ Check[val /. replacementRules, val]],
            expr,
            10
        ];
        allFlowVars = Join[Lookup[Eqs, "js", {}], Lookup[Eqs, "jts", {}]];
        resolved = Association @ Cases[
            Map[
                Function[{var},
                    Module[{val},
                        val = Quiet @ Check[N @ resolve[var], Missing["NotAvailable"]];
                        If[NumericQ[val],
                            var -> If[Abs[val] <= tol, 0., N[val]],
                            Nothing
                        ]
                    ]
                ],
                allFlowVars
            ],
            _Rule
        ];
        If[
            !AssociationQ[resolved] ||
            !And @@ (KeyExistsQ[resolved, #] & /@ allFlowVars),
            CriticalJFirstFailure["FlowReconstructionFailed"],
            resolved
        ]
    ];

SolveCriticalJFirstUtilities[
    Eqs_Association,
    flowAssoc_Association,
    uVars_List,
    tol_?NumericQ
] :=
    Module[{costpluscurrents, eqGeneralCritical, eqValueAuxiliaryEdges, ruleEntryValues,
      ruleExitValues, ineqSwitchingByVertex, switchingEqualities, uEqualities,
      uSystem, uCandidate, residual, utilityAssoc},
        If[uVars === {},
            Return[CriticalJFirstFailure["MissingUVars"], Module]
        ];
        costpluscurrents = Lookup[Eqs, "costpluscurrents", <||>];
        eqGeneralCritical =
            Lookup[Eqs, "EqGeneral", True] /.
                If[AssociationQ[costpluscurrents],
                    Expand /@ (costpluscurrents /. flowAssoc),
                    {}
                ];
        eqValueAuxiliaryEdges = Lookup[Eqs, "EqValueAuxiliaryEdges", True];
        ruleEntryValues = Lookup[Eqs, "RuleEntryValues", <||>];
        ruleExitValues = Lookup[Eqs, "RuleExitValues", <||>];
        ineqSwitchingByVertex = Lookup[Eqs, "IneqSwitchingByVertex", True];
        switchingEqualities = DeleteDuplicatesBy[
            Cases[
                CriticalConjunctionTerms[ineqSwitchingByVertex],
                LessEqual[lhs_, rhs_] :> (lhs == rhs)
            ],
            Sort @* (List @@ # &)
        ];
        uEqualities = DeleteDuplicates @ Select[
            Cases[
                DeleteCases[
                    Flatten[{
                        CriticalConjunctionTerms[eqGeneralCritical],
                        CriticalConjunctionTerms[eqValueAuxiliaryEdges],
                        switchingEqualities,
                        KeyValueMap[Equal, ruleEntryValues],
                        KeyValueMap[Equal, ruleExitValues]
                    }] /. flowAssoc,
                    True
                ],
                _Equal
            ],
            !FreeQ[#, Alternatives @@ uVars] &
        ];
        If[uEqualities === {},
            Return[CriticalJFirstFailure["URecoveryFailed"], Module]
        ];
        uSystem = BuildLinearSystemFromEqualities[uEqualities, uVars];
        If[!AssociationQ[uSystem],
            Return[CriticalJFirstFailure["URecoveryFailed"], Module]
        ];
        uCandidate = LinearSolveCandidate[uSystem["A"], uSystem["b"]];
        If[!VectorQ[uCandidate, NumericQ] || Length[uCandidate] =!= Length[uVars],
            Return[CriticalJFirstFailure["URecoveryFailed"], Module]
        ];
        residual = Quiet @ Check[N @ Norm[uSystem["A"] . uCandidate - uSystem["b"], Infinity], Infinity];
        If[!NumericQ[residual] || residual > tol,
            Return[CriticalJFirstFailure["URecoveryFailed"], Module]
        ];
        utilityAssoc = AssociationThread[uVars, uCandidate];
        Association @ KeyValueMap[
            #1 -> If[Abs[#2] <= tol, 0., N[#2]] &,
            utilityAssoc
        ]
    ];

SolveCriticalJFirstBackend[Eqs_Association, numericState_Association] :=
    Module[{decoupling, topo, massConstraints, massVars, aEq, bEqRaw, bEq, tol = $CriticalSolverTolerance,
      nVars, massAssoc = <||>, xCandidate, candidateResidual, flowMin, rank,
      objectiveExpr, qpVars, qpRules, objectiveQP, constraintsQP,
      optimizerOrder, optimizer, rawSolution = $Failed, massVec,
      flowAssoc, uVars, utilityAssoc, solvedAssoc, expectedVars},
        decoupling = Lookup[numericState, "CriticalDecoupling", Missing["NotAvailable"]];
        If[!AssociationQ[decoupling],
            Return[CriticalJFirstFailure["MissingDecoupling"], Module]
        ];
        topo = Lookup[decoupling, "TopologicalOrder", <||>];
        If[!TrueQ[Lookup[topo, "IsDAG", False]],
            Return[CriticalJFirstFailure["NonDAG"], Module]
        ];
        massConstraints = Lookup[decoupling, "MassConstraints", <||>];
        massVars = Lookup[massConstraints, "Variables", {}];
        aEq = Lookup[massConstraints, "Aeq", Missing["NotAvailable"]];
        bEqRaw = Lookup[massConstraints, "beq", Missing["NotAvailable"]];
        If[Head[aEq] =!= SparseArray,
            aEq = SparseArray[aEq]
        ];
        bEq = Developer`ToPackedArray @ N @ Which[
            Head[bEqRaw] === SparseArray, Normal[bEqRaw],
            VectorQ[bEqRaw, NumericQ], bEqRaw,
            True, {}
        ];
        nVars = Length[massVars];
        If[
            nVars == 0 || Head[aEq] =!= SparseArray || !VectorQ[bEq, NumericQ] ||
            Length[bEq] =!= Dimensions[aEq][[1]] || Dimensions[aEq][[2]] =!= nVars,
            Return[CriticalJFirstFailure["MissingMassConstraints"], Module]
        ];

        rank = Quiet @ Check[MatrixRank[Normal[aEq]], -1];
        If[IntegerQ[rank] && rank === nVars,
            xCandidate = LinearSolveCandidate[aEq, bEq];
            If[VectorQ[xCandidate, NumericQ] && Length[xCandidate] === nVars,
                candidateResidual = Quiet @ Check[N @ Norm[aEq . xCandidate - bEq, Infinity], Infinity];
                flowMin = Min[xCandidate];
                If[
                    NumericQ[candidateResidual] && candidateResidual <= tol &&
                    NumericQ[flowMin] && flowMin >= -tol,
                    massAssoc = AssociationThread[massVars, xCandidate]
                ]
            ]
        ];

        If[!AssociationQ[massAssoc] || massAssoc === <||>,
            objectiveExpr = BuildCriticalJFirstObjectiveExpression[Eqs, massVars];
            qpVars = Table[Unique["qv"], {nVars}];
            qpRules = Thread[massVars -> qpVars];
            constraintsQP = Join[
                Thread[(aEq . qpVars) == bEq],
                Thread[qpVars >= 0]
            ];
            If[!MatchQ[objectiveExpr, Failure[_, _Association]],
                objectiveQP = objectiveExpr /. qpRules;
                optimizerOrder = DeleteDuplicates @ {
                    If[
                        NameQ["System`QuadraticOptimization"] &&
                        TrueQ @ Quiet @ Check[UseQuadraticCriticalBackendQ[Eqs], False],
                        QuadraticOptimization,
                        ConvexOptimization
                    ],
                    ConvexOptimization
                };
                Do[
                    rawSolution = Quiet @ Check[
                        optimizer[
                            objectiveQP,
                            constraintsQP,
                            qpVars
                        ],
                        $Failed
                    ];
                    If[MatchQ[rawSolution, {_Rule ..}],
                        Break[]
                    ],
                    {optimizer, optimizerOrder}
                ];
            ];
            If[!MatchQ[rawSolution, {_Rule ..}],
                rawSolution = Quiet @ Check[
                    LinearOptimization[
                        Total[qpVars],
                        constraintsQP,
                        qpVars
                    ],
                    $Failed
                ];
            ];
            If[!MatchQ[rawSolution, {_Rule ..}],
                Return[CriticalJFirstFailure["FlowInfeasible"], Module]
            ];
            massAssoc = AssociationThread[massVars, N @ (qpVars /. rawSolution)];
        ];

        massVec = Lookup[massAssoc, massVars, Missing["NotAvailable"]];
        If[!VectorQ[massVec, NumericQ] || Length[massVec] =!= nVars,
            Return[CriticalJFirstFailure["FlowInfeasible"], Module]
        ];
        candidateResidual = Quiet @ Check[N @ Norm[aEq . massVec - bEq, Infinity], Infinity];
        flowMin = Min[massVec];
        If[
            !NumericQ[candidateResidual] || candidateResidual > tol ||
            !NumericQ[flowMin] || flowMin < -tol,
            Return[CriticalJFirstFailure["FlowInfeasible"], Module]
        ];

        flowAssoc = RecoverCriticalFlowAssociation[Eqs, massAssoc, tol];
        If[MatchQ[flowAssoc, Failure[_, _Association]],
            Return[flowAssoc, Module]
        ];
        uVars = Lookup[decoupling, "UVars", Lookup[numericState, "UVars", {}]];
        utilityAssoc = SolveCriticalJFirstUtilities[Eqs, flowAssoc, uVars, tol];
        If[MatchQ[utilityAssoc, Failure[_, _Association]],
            Return[utilityAssoc, Module]
        ];
        solvedAssoc = Join[utilityAssoc, flowAssoc];
        expectedVars = Join[
            Lookup[Eqs, "us", {}],
            Lookup[Eqs, "js", {}],
            Lookup[Eqs, "jts", {}]
        ];
        If[!And @@ (KeyExistsQ[solvedAssoc, #] & /@ expectedVars),
            Return[CriticalJFirstFailure["IncompleteCandidate"], Module]
        ];
        Join[
            KeyTake[solvedAssoc, Lookup[Eqs, "us", {}]],
            KeyTake[solvedAssoc, Lookup[Eqs, "js", {}]],
            KeyTake[solvedAssoc, Lookup[Eqs, "jts", {}]]
        ]
    ];

SolveCriticalNumericBackendWithTelemetry[Eqs_Association] :=
    Module[{numericState, jFirstCandidate, jFirstReason, coupledCandidate,
      jFirstStatus, jFirstValidQ, forcedRejectQ},
        numericState = Lookup[Eqs, "NumericState", <||>];
        forcedRejectQ = TrueQ[Lookup[Eqs, "ForceNumericBackendValidationFailure", False]];
        jFirstCandidate = SolveCriticalJFirstBackend[Eqs, numericState];
        If[AssociationQ[jFirstCandidate],
            jFirstStatus = CheckFlowFeasibility[jFirstCandidate];
            jFirstValidQ =
                TrueQ[
                    !forcedRejectQ &&
                    jFirstStatus === "Feasible" &&
                    IsCriticalSolution[
                        Join[Eqs, <|"AssoCritical" -> jFirstCandidate, "Solution" -> jFirstCandidate|>],
                        "Tolerance" -> $CriticalSolverTolerance,
                        "Verbose" -> False
                    ]
                ];
            If[jFirstValidQ,
                Return[
                    <|
                        "Candidate" -> jFirstCandidate,
                        "Strategy" -> "JFirst",
                        "JFirstBackendUsed" -> True,
                        "JFirstBackendFallbackReason" -> None
                    |>,
                    Module
                ]
            ];
            jFirstReason = Which[
                forcedRejectQ, "ForcedValidationFailure",
                jFirstStatus =!= "Feasible", "FlowFeasibilityFailed",
                True, "CriticalValidationFailed"
            ],
            jFirstReason = CriticalJFirstFailureReason[jFirstCandidate]
        ];
        coupledCandidate = SolveCriticalNumericBackend[Eqs];
        <|
            "Candidate" -> coupledCandidate,
            "Strategy" -> "Coupled",
            "JFirstBackendUsed" -> False,
            "JFirstBackendFallbackReason" -> jFirstReason
        |>
    ];

SolveCriticalNumericBackend[Eqs_Association] :=
    Module[{constraints, equalities, flowVars, allVars, fullUs, fullJs, fullJts,
        reduction, subRules, reducedUs, reducedVars,
        linearSystem, aMat, bVec, xCandidate, varIndex, flowIndices,
        linResidual, flowMin, tol = $CriticalSolverTolerance, assoc, lpResult, lpVals,
        solvedAssoc, fullAssoc,
        nVars, cVec, aIneq, bIneq, aEq, bEq},
        constraints = BuildCriticalLinearConstraints[Eqs];
        If[constraints === $Failed || !AssociationQ[constraints],
            Return[$Failed, Module]
        ];
        equalities = Lookup[constraints, "Equalities", {}];
        flowVars = Lookup[constraints, "FlowVariables", {}];
        allVars = Lookup[constraints, "AllVariables", {}];
        fullUs = Lookup[Eqs, "us", {}];
        fullJs = Lookup[Eqs, "js", {}];
        fullJts = Lookup[Eqs, "jts", {}];
        If[equalities === {} || allVars === {},
            Return[$Failed, Module]
        ];

        (* Apply u-variable equivalence reduction *)
        reduction = BuildUEquivalenceReduction[Eqs];
        subRules = Lookup[reduction, "SubstitutionRules", {}];
        reducedUs = Lookup[reduction, "ReducedUVars", fullUs];
        If[subRules =!= {},
            equalities = DeleteDuplicates @ Cases[
                DeleteCases[equalities /. subRules, True],
                _Equal
            ];
            (* Remove tautologies like x == x that arise from substitution *)
            equalities = Select[equalities, #[[1]] =!= #[[2]] &];
            reducedVars = Join[reducedUs, fullJs, fullJts];
            ,
            reducedVars = allVars
        ];
        If[equalities === {} || reducedVars === {},
            Return[$Failed, Module]
        ];

        varIndex = AssociationThread[reducedVars, Range[Length[reducedVars]]];
        flowIndices = Lookup[varIndex, flowVars, Missing["NotAvailable"]];
        If[MemberQ[flowIndices, _Missing] || flowIndices === {},
            Return[$Failed, Module]
        ];

        linearSystem = BuildLinearSystemFromEqualities[equalities, reducedVars];
        If[AssociationQ[linearSystem],
            aMat = linearSystem["A"];
            bVec = linearSystem["b"];
            xCandidate = LinearSolveCandidate[aMat, bVec];
            If[VectorQ[xCandidate, NumericQ] && Length[xCandidate] === Length[reducedVars],
                linResidual = Quiet @ Check[N @ Norm[aMat . xCandidate - bVec, Infinity], Infinity];
                flowMin = Min[xCandidate[[flowIndices]]];
                If[NumericQ[linResidual] && linResidual <= tol && NumericQ[flowMin] && flowMin >= -tol,
                    assoc = AssociationThread[reducedVars, xCandidate];
                    (* Lift eliminated u-vars back from their representatives *)
                    fullAssoc = Join[assoc, Association[subRules /. assoc]];
                    solvedAssoc = Association @ KeyValueMap[
                        #1 -> If[Abs[#2] <= tol, 0., N[#2]] &,
                        fullAssoc
                    ];
                    Return[
                        Join[
                            KeyTake[solvedAssoc, fullUs],
                            KeyTake[solvedAssoc, fullJs],
                            KeyTake[solvedAssoc, fullJts]
                        ],
                        Module
                    ]
                ]
            ]
            ,
            Return[$Failed, Module]
        ];

        nVars = Length[reducedVars];
        (* Infinitesimal flow-minimizing objective to break counter-flow
           degeneracy without distorting the equilibrium.  The equilibrium
           is fully pinned by the equality constraints; this epsilon
           tiebreaker only selects the non-cycling vertex of the polytope. *)
        cVec = Developer`ToPackedArray @ N @ ReplacePart[
            ConstantArray[0., nVars],
            Thread[flowIndices -> 10.^-12]
        ];
        aEq = If[Head[aMat] === SparseArray, aMat, SparseArray[aMat]];
        bEq = Developer`ToPackedArray @ N @ (-bVec);
        aIneq = SparseArray[
            Thread[
                Transpose[{Range[Length[flowIndices]], flowIndices}] -> 1.
            ],
            {Length[flowIndices], nVars}
        ];
        bIneq = Developer`ToPackedArray @ ConstantArray[0., Length[flowIndices]];
        lpResult = Quiet @ Check[
            LinearOptimization[
                cVec,
                {aIneq, bIneq},
                {aEq, bEq},
                Tolerance -> tol
            ],
            $Failed
        ];
        lpVals = If[
            VectorQ[lpResult, NumericQ] && Length[lpResult] === nVars,
            Developer`ToPackedArray @ N @ lpResult,
            $Failed
        ];
        If[lpVals === $Failed,
            lpResult = Quiet @ Check[
                LinearOptimization[
                    cVec,
                    {aIneq, bIneq},
                    {aEq, bEq},
                    Method -> "Simplex",
                    Tolerance -> tol
                ],
                $Failed
            ];
            lpVals = If[
                VectorQ[lpResult, NumericQ] && Length[lpResult] === nVars,
                Developer`ToPackedArray @ N @ lpResult,
                $Failed
            ];
        ];
        If[
            lpVals === $Failed || !VectorQ[lpVals, NumericQ] || Length[lpVals] =!= nVars,
            Return[$Failed, Module]
        ];
        assoc = AssociationThread[reducedVars, Developer`ToPackedArray @ N @ lpVals];
        (* Lift eliminated u-vars back from their representatives *)
        fullAssoc = Join[assoc, Association[subRules /. assoc]];
        solvedAssoc = Association @ KeyValueMap[
            #1 -> If[Abs[#2] <= tol, 0., N[#2]] &,
            fullAssoc
        ];
        Join[
            KeyTake[solvedAssoc, fullUs],
            KeyTake[solvedAssoc, fullJs],
            KeyTake[solvedAssoc, fullJts]
        ]
    ];


BuildCriticalResult[PreEqs_Association, resultKind_, status_, message_,
    candidate_, telemetry_Association, unresolvedEqs_: True] :=
    Module[{solution, comparisonData},
        solution = If[AssociationQ[candidate], candidate, Missing["NotAvailable"]];
        comparisonData = BuildSolverComparisonData[PreEqs, solution];
        Join[PreEqs, MakeSolverResult["CriticalCongestion", resultKind,
            status, message, solution,
            Join[comparisonData,
                <|"AssoCritical" -> candidate, "Status" -> status,
                  "UnresolvedEquations" -> unresolvedEqs|>,
                telemetry]
        ]]
    ];

(* EnsurePreprocessed: lazily run MFGPreprocessing if InitRules are not already present. *)
EnsurePreprocessed[Eqs_Association] :=
    Module[{time, PreEqs},
        If[KeyExistsQ[Eqs, "InitRules"],
            Eqs,
            {time, PreEqs} = AbsoluteTiming @ MFGPreprocessing[Eqs];
            MFGPrint["Preprocessing (for downstream solvers) took ", time, " seconds."];
            PreEqs
        ]
    ];

Options[CriticalCongestionSolver] = {
    "ValidateSolution" -> True,
    "ValidationTolerance" -> $CriticalSolverTolerance,
    "ValidationVerbose" -> False,
    "SymbolicTimeLimit" -> 120.,
    "ExactMode" -> False
};

CriticalCongestionSolver::noassoc = "Expected an Association for Eqs, but received `1`. Ensure input is pre-processed via DataToEquations.";

CriticalCongestionSolver[$Failed, ___] :=
    $Failed

CriticalCongestionSolver[Eqs : Except[_Association | $Failed], OptionsPattern[]] :=
    (Message[CriticalCongestionSolver::noassoc, Head[Eqs]]; $Failed)

CriticalCongestionSolver[Eqs_Association, OptionsPattern[]] :=
    Module[{PreEqs, js, AssoCritical, time, status,
         resultKind, message,
         validateQ, validationTolerance, validationVerboseQ, validationFailedQ,
         numericBackendRequestedQ, numericBackendEligibleQ,
         numericBackendUsedQ, numericBackendFallbackReason, numericBackendSolveTime,
         numericBackendStrategy, jFirstBackendUsedQ, jFirstBackendFallbackReason,
         numericCandidate, numericOutcome, forcedRejectQ, numericBackendTimeLimit,
         symbolicTimeLimit, exactModeQ,
         telemetry},
        ClearSolveCache[];
        exactModeQ = TrueQ[OptionValue["ExactMode"]];
        validateQ = TrueQ[OptionValue["ValidateSolution"]];
        validationTolerance = OptionValue["ValidationTolerance"];
        validationVerboseQ = TrueQ[OptionValue["ValidationVerbose"]];
        validationFailedQ = False;
        numericBackendUsedQ = False;
        numericBackendSolveTime = Missing["NotAvailable"];
        numericBackendStrategy = "Symbolic";
        jFirstBackendUsedQ = False;
        jFirstBackendFallbackReason = None;
        numericBackendRequestedQ = !exactModeQ && CriticalNumericBackendRequestedQ[Eqs];
        numericBackendEligibleQ = !exactModeQ && CriticalNumericBackendEligibleQ[Eqs];
        numericBackendFallbackReason = Which[
            !numericBackendRequestedQ, "DisabledByGuard",
            !numericBackendEligibleQ, "IneligibleInput",
            True, "NumericSolveFailed"
        ];
        forcedRejectQ = TrueQ[Lookup[Eqs, "ForceNumericBackendValidationFailure", False]];
        numericBackendTimeLimit = Lookup[Eqs, "CriticalNumericBackendTimeLimit", 30.];
        numericBackendTimeLimit = Which[
            numericBackendTimeLimit === Infinity || numericBackendTimeLimit === DirectedInfinity[1],
                Infinity,
            NumericQ[numericBackendTimeLimit] && numericBackendTimeLimit > 0,
                N[numericBackendTimeLimit],
            True,
                30.
        ];
        symbolicTimeLimit = OptionValue["SymbolicTimeLimit"];
        symbolicTimeLimit = Which[
            symbolicTimeLimit === Infinity || symbolicTimeLimit === DirectedInfinity[1],
                Infinity,
            NumericQ[symbolicTimeLimit] && symbolicTimeLimit > 0,
                N[symbolicTimeLimit],
            True,
                120.
        ];

        (* Helper: build telemetry association from current state *)
        telemetry := <|
            "NumericBackendUsed" -> numericBackendUsedQ,
            "NumericBackendFallbackReason" -> numericBackendFallbackReason,
            "NumericBackendSolveTime" -> numericBackendSolveTime,
            "NumericBackendStrategy" -> numericBackendStrategy,
            "JFirstBackendUsed" -> jFirstBackendUsedQ,
            "JFirstBackendFallbackReason" -> jFirstBackendFallbackReason
        |>;

        (* Try the numeric backend first when explicitly forced or globally enabled. *)
        If[numericBackendRequestedQ && numericBackendEligibleQ,
            MFGPrint["Attempting critical numeric backend..."];
            {numericBackendSolveTime, numericOutcome} =
                AbsoluteTiming[
                    TimeConstrained[
                        SolveCriticalNumericBackendWithTelemetry[Eqs],
                        numericBackendTimeLimit,
                        $TimedOut
                    ]
                ];
            If[numericOutcome === $TimedOut,
                numericBackendFallbackReason = "NumericBackendTimeout";
                numericBackendStrategy = "Symbolic";
                jFirstBackendUsedQ = False;
                jFirstBackendFallbackReason = "NumericBackendTimeout";
                MFGPrint[
                    "Numeric backend timed out after ",
                    numericBackendTimeLimit,
                    " seconds; falling back to symbolic path."
                ];
                numericCandidate = $Failed;
                ,
                numericCandidate = Lookup[numericOutcome, "Candidate", $Failed];
                numericBackendStrategy = Lookup[numericOutcome, "Strategy", "Coupled"];
                jFirstBackendUsedQ = TrueQ[Lookup[numericOutcome, "JFirstBackendUsed", False]];
                jFirstBackendFallbackReason = Lookup[
                    numericOutcome,
                    "JFirstBackendFallbackReason",
                    "NotAttempted"
                ];
            ];
            If[AssociationQ[numericCandidate] && numericCandidate =!= $Failed,
                status = CheckFlowFeasibility[numericCandidate];
                If[status === "Feasible",
                    numericBackendUsedQ = True;
                    numericBackendFallbackReason = None;
                    PreEqs = EnsurePreprocessed[Eqs];
                    If[validateQ,
                        Module[{valReport},
                            valReport = If[forcedRejectQ,
                                <|"Valid" -> False, "Reason" -> "ForcedValidationFailure",
                                  "IndeterminateResiduals" -> <||>|>,
                                IsCriticalSolution[
                                    Join[PreEqs, <|"AssoCritical" -> numericCandidate, "Solution" -> numericCandidate|>],
                                    "Tolerance" -> validationTolerance,
                                    "Verbose" -> validationVerboseQ,
                                    "ReturnReport" -> True
                                ]
                            ];
                            Which[
                                valReport["Reason"] === "IndeterminateConstraints",
                                    Return[BuildCriticalResult[PreEqs, "Success", status, None,
                                        numericCandidate, telemetry,
                                        Lookup[valReport, "IndeterminateResiduals", <||>]], Module]
                                ,
                                !TrueQ[valReport["Valid"]],
                                    numericBackendUsedQ = False;
                                    numericBackendStrategy = "Symbolic";
                                    If[
                                        TrueQ[jFirstBackendUsedQ] &&
                                        (jFirstBackendFallbackReason === None || jFirstBackendFallbackReason === "NotAttempted"),
                                        jFirstBackendFallbackReason =
                                            If[forcedRejectQ, "ForcedValidationFailure", "CriticalValidationFailed"]
                                    ];
                                    jFirstBackendUsedQ = False;
                                    numericBackendFallbackReason =
                                        If[forcedRejectQ, "ForcedValidationFailure", "CriticalValidationFailed"]
                                ,
                                True,
                                    Return[BuildCriticalResult[PreEqs, "Success", status, None, numericCandidate, telemetry], Module]
                            ]
                        ]
                        ,
                        Return[BuildCriticalResult[PreEqs, "Success", status, None, numericCandidate, telemetry], Module]
                    ],
                    numericBackendFallbackReason = "FlowFeasibilityFailed"
                ],
                If[numericBackendFallbackReason =!= "NumericBackendTimeout",
                    numericBackendFallbackReason = "NumericSolveFailed"
                ]
            ]
        ];

        (* Try direct solver for zero-switching-cost, fully numeric networks *)
        If[!exactModeQ &&
            AllTrue[
                Values[Lookup[Eqs, "SwitchingCosts", <|_ -> 1|>]],
                NumericQ[#] && PossibleZeroQ[#] &
            ] &&
            Lookup[Eqs, "auxiliaryGraph", None] =!= None &&
            AllTrue[Last /@ Lookup[Eqs, "Entrance Vertices and Flows", {}], NumericQ] &&
            AllTrue[Last /@ Lookup[Eqs, "Exit Vertices and Terminal Costs", {}], NumericQ] &&
            GraphDistanceHeuristicSafeQ[Eqs],
            MFGPrint["Attempting direct critical solver (zero switching costs)..."];
            {time, AssoCritical} = AbsoluteTiming[DirectCriticalSolver[Eqs]];
            If[AssociationQ[AssoCritical] && AssoCritical =!= $Failed,
                MFGPrint["Direct critical solver completed in ", time, " seconds."];
                status = CheckFlowFeasibility[AssoCritical];
                If[status === "Feasible",
                    PreEqs = EnsurePreprocessed[Eqs];
                    If[validateQ,
                        Module[{valReport},
                            valReport = IsCriticalSolution[
                                Join[PreEqs, <|"AssoCritical" -> AssoCritical, "Solution" -> AssoCritical|>],
                                "Tolerance" -> validationTolerance,
                                "Verbose"   -> validationVerboseQ,
                                "ReturnReport" -> True
                            ];
                            Which[
                                valReport["Reason"] === "IndeterminateConstraints",
                                    Return[BuildCriticalResult[PreEqs, "Success", status, None,
                                        AssoCritical, telemetry,
                                        Lookup[valReport, "IndeterminateResiduals", <||>]], Module]
                                ,
                                !TrueQ[valReport["Valid"]],
                                    MFGPrint["Direct solver validation failed, falling back to symbolic solver."]
                                ,
                                True,
                                    Return[BuildCriticalResult[PreEqs, "Success", status, None, AssoCritical, telemetry], Module]
                            ]
                        ]
                        ,
                        Return[BuildCriticalResult[PreEqs, "Success", status, None, AssoCritical, telemetry], Module]
                    ],
                    MFGPrint["Direct solver produced infeasible result, falling back to symbolic solver."]
                ],
                MFGPrint["Direct solver failed, falling back to symbolic solver."]
            ]
        ];
        (* Standard symbolic solver path *)
        PreEqs = EnsurePreprocessed[If[AssociationQ[PreEqs], PreEqs, Eqs]];
        (* Oracle Bridge: if the numeric backend produced a candidate that failed validation,
           use it to prune the symbolic Or-branches before the full symbolic solve. *)
        If[!exactModeQ &&
            AssociationQ[numericCandidate] &&
            ListQ[Lookup[PreEqs, "NewSystem", $Failed]] &&
            Length[Lookup[PreEqs, "NewSystem", {}]] === 3,
            Module[{prunedSystem, originalBranchCount, prunedBranchCount},
                originalBranchCount = Length[Lookup[PreEqs, "NewSystem"][[3]]];
                prunedSystem = BuildPrunedSystem[
                    Lookup[PreEqs, "NewSystem"],
                    numericCandidate
                ];
                prunedBranchCount = Length[prunedSystem[[3]]];
                If[prunedBranchCount < originalBranchCount,
                    PreEqs = Append[PreEqs, "NewSystem" -> prunedSystem];
                    MFGPrint["Oracle Bridge: pruned Or-branches from ", originalBranchCount,
                        " to ", prunedBranchCount]
                ]
            ]
        ];
        js = Lookup[PreEqs, "js", $Failed];
        Module[{symbolicResult},
            symbolicResult = TimeConstrained[
                Module[{localAssoCritical, localStatus, localValidationFailedQ = False,
                        localResultKind, localMessage, mfgssResult, unresolvedUConstraints, valReport,
                        symbolicTelemetry},
                    mfgssResult = MFGSystemSolver[PreEqs][AssociationThread[js, ConstantArray[0, Length[js]]]];
                    localAssoCritical = mfgssResult["Solution"];
                    unresolvedUConstraints = mfgssResult["UnresolvedConstraints"];
                    symbolicTelemetry = Append[telemetry, "SymbolicRegion" -> unresolvedUConstraints];
                    localStatus = CheckFlowFeasibility[localAssoCritical];
                    If[validateQ && AssociationQ[localAssoCritical] && localStatus === "Feasible",
                        valReport = IsCriticalSolution[
                            Join[PreEqs, <|"AssoCritical" -> localAssoCritical, "Solution" -> localAssoCritical|>],
                            "Tolerance" -> validationTolerance,
                            "Verbose"   -> validationVerboseQ,
                            "ReturnReport" -> True
                        ];
                        Which[
                            valReport["Reason"] === "IndeterminateConstraints",
                                Return[BuildCriticalResult[PreEqs, "Success", localStatus, None,
                                    localAssoCritical, symbolicTelemetry,
                                    If[unresolvedUConstraints =!= None,
                                        Join[Lookup[valReport, "IndeterminateResiduals", <||>],
                                             <|"FlowSolveResidual" -> unresolvedUConstraints|>],
                                        Lookup[valReport, "IndeterminateResiduals", <||>]
                                    ]], Module]
                            ,
                            !TrueQ[valReport["Valid"]],
                                localValidationFailedQ = True;
                                localAssoCritical = Null;
                                localStatus = "Infeasible"
                            ,
                            True, Null (* valid — fall through *)
                        ]
                    ];
                    localResultKind = If[localAssoCritical === Null, "Failure", "Success"];
                    localMessage = If[localAssoCritical === Null,
                        If[localValidationFailedQ, "InvalidCriticalSolution", "NoSolution"], None];
                    BuildCriticalResult[PreEqs, localResultKind, localStatus, localMessage, localAssoCritical, symbolicTelemetry]
                ],
                symbolicTimeLimit,
                $TimedOut
            ];
            If[symbolicResult === $TimedOut,
                BuildCriticalResult[PreEqs, "Failure", Missing["NotAvailable"],
                    "SymbolicTimeout", Null,
                    Join[telemetry, <|"SymbolicSolverTimedOut" -> True,
                        "SymbolicSolverTimeLimit" -> symbolicTimeLimit|>]]
                ,
                symbolicResult
            ]
        ]
    ];

Options[IsCriticalSolution] = {
    "Tolerance" -> $CriticalSolverTolerance,
    "Verbose" -> False,
    "ReturnReport" -> False
};

NumericRelationSatisfiedQ[expr_, tol_?NumericQ] :=
    Module[{delta, pair, held},
        held = HoldComplete[expr];
        Which[
            MatchQ[held, HoldComplete[Equal[_, _]] | HoldComplete[Inactive[Equal][_, _]]],
                pair = Replace[held, {
                    HoldComplete[Equal[l_, r_]] :> {l, r},
                    HoldComplete[Inactive[Equal][l_, r_]] :> {l, r}
                }];
                delta = Quiet @ Check[N[pair[[1]] - pair[[2]]], $Failed];
                Which[delta === $Failed, False, NumericQ[delta], Abs[delta] <= tol, True, Indeterminate]
            ,
            MatchQ[held, HoldComplete[LessEqual[_, _]] | HoldComplete[Inactive[LessEqual][_, _]]],
                pair = Replace[held, {
                    HoldComplete[LessEqual[l_, r_]] :> {l, r},
                    HoldComplete[Inactive[LessEqual][l_, r_]] :> {l, r}
                }];
                delta = Quiet @ Check[N[pair[[1]] - pair[[2]]], $Failed];
                Which[delta === $Failed, False, NumericQ[delta], delta <= tol, True, Indeterminate]
            ,
            MatchQ[held, HoldComplete[GreaterEqual[_, _]] | HoldComplete[Inactive[GreaterEqual][_, _]]],
                pair = Replace[held, {
                    HoldComplete[GreaterEqual[l_, r_]] :> {l, r},
                    HoldComplete[Inactive[GreaterEqual][l_, r_]] :> {l, r}
                }];
                delta = Quiet @ Check[N[pair[[2]] - pair[[1]]], $Failed];
                Which[delta === $Failed, False, NumericQ[delta], delta <= tol, True, Indeterminate]
            ,
            MatchQ[held, HoldComplete[Less[_, _]] | HoldComplete[Inactive[Less][_, _]]],
                pair = Replace[held, {
                    HoldComplete[Less[l_, r_]] :> {l, r},
                    HoldComplete[Inactive[Less][l_, r_]] :> {l, r}
                }];
                delta = Quiet @ Check[N[pair[[1]] - pair[[2]]], $Failed];
                Which[delta === $Failed, False, NumericQ[delta], delta < -tol, True, Indeterminate]
            ,
            MatchQ[held, HoldComplete[Greater[_, _]] | HoldComplete[Inactive[Greater][_, _]]],
                pair = Replace[held, {
                    HoldComplete[Greater[l_, r_]] :> {l, r},
                    HoldComplete[Inactive[Greater][l_, r_]] :> {l, r}
                }];
                delta = Quiet @ Check[N[pair[[2]] - pair[[1]]], $Failed];
                Which[delta === $Failed, False, NumericQ[delta], delta < -tol, True, Indeterminate]
            ,
            MatchQ[held, HoldComplete[Unequal[_, _]] | HoldComplete[Inactive[Unequal][_, _]]],
                pair = Replace[held, {
                    HoldComplete[Unequal[l_, r_]] :> {l, r},
                    HoldComplete[Inactive[Unequal][l_, r_]] :> {l, r}
                }];
                delta = Quiet @ Check[N[pair[[1]] - pair[[2]]], $Failed];
                Which[delta === $Failed, False, NumericQ[delta], Abs[delta] > tol, True, Indeterminate]
            ,
            True,
                False
        ]
    ];

LogicalSatisfiedQ[expr_, rules_Association, tol_?NumericQ] :=
    Module[{heldExpr, evaluated},
        Which[
            expr === True, True,
            expr === False, False,
            True,
                heldExpr =
                    expr /. HoldPattern[(h:Equal | LessEqual | GreaterEqual | Less | Greater | Unequal)[lhs_, rhs_]] :>
                        Inactive[h][lhs, rhs];
                evaluated = Quiet @ Check[
                    heldExpr /. rules /. {
                        HoldPattern[Inactive[Equal][lhs_, rhs_]] :>
                            NumericRelationSatisfiedQ[Inactive[Equal][lhs, rhs], tol],
                        HoldPattern[Inactive[LessEqual][lhs_, rhs_]] :>
                            NumericRelationSatisfiedQ[Inactive[LessEqual][lhs, rhs], tol],
                        HoldPattern[Inactive[GreaterEqual][lhs_, rhs_]] :>
                            NumericRelationSatisfiedQ[Inactive[GreaterEqual][lhs, rhs], tol],
                        HoldPattern[Inactive[Less][lhs_, rhs_]] :>
                            NumericRelationSatisfiedQ[Inactive[Less][lhs, rhs], tol],
                        HoldPattern[Inactive[Greater][lhs_, rhs_]] :>
                            NumericRelationSatisfiedQ[Inactive[Greater][lhs, rhs], tol],
                        HoldPattern[Inactive[Unequal][lhs_, rhs_]] :>
                            NumericRelationSatisfiedQ[Inactive[Unequal][lhs, rhs], tol]
                    },
                    $Failed
                ];
                If[evaluated === $Failed,
                    False,
                    Module[{simplified = Simplify[evaluated]},
                        Which[
                            TrueQ[simplified],    True,
                            simplified === False,  False,
                            True,                  Indeterminate
                        ]
                    ]
                ]
        ]
    ];

RulesConjunction[rules_] :=
    Module[{pairs},
        pairs = Which[
            AssociationQ[rules], Normal[rules],
            ListQ[rules] && VectorQ[rules, MatchQ[#, _Rule]&], rules,
            True, {}
        ];
        If[pairs === {},
            True,
            And @@ (Equal @@@ (List @@@ pairs))
        ]
    ];

CriticalSolutionAssociation[Eqs_Association] :=
    Module[{candidate},
        candidate = Lookup[Eqs, "AssoCritical", Missing["NotAvailable"]];
        If[AssociationQ[candidate],
            candidate,
            candidate = Lookup[Eqs, "Solution", Missing["NotAvailable"]];
            If[AssociationQ[candidate], candidate, Missing["NotAvailable"]]
        ]
    ];

CriticalFlowApproximation[Eqs_Association, solution_Association] :=
    KeyTake[solution, Lookup[Eqs, "js", {}]];

IsCriticalSolution[Eqs_Association, OptionsPattern[]] :=
    Module[{tol, verboseQ, returnReportQ, solution, blocks, blockResults,
         residual, residualPass, overall, report, flowApprox,
         costpluscurrents, eqGeneralCritical, eqResidualPairs, eqResidualDeltas,
         requiredKeys, missingKeys,
         concretelyFailed, indeterminateBlocks, hasIndeterminate, indeterminateResiduals},
        tol = OptionValue["Tolerance"];
        verboseQ = TrueQ[OptionValue["Verbose"]];
        returnReportQ = TrueQ[OptionValue["ReturnReport"]];

        requiredKeys = {
            "EqEntryIn",
            "EqValueAuxiliaryEdges",
            "AltOptCond",
            "EqBalanceSplittingFlows",
            "AltFlows",
            "AltTransitionFlows",
            "IneqJs",
            "IneqJts",
            "IneqSwitchingByVertex",
            "RuleEntryOut",
            "RuleExitFlowsIn",
            "RuleEntryValues",
            "RuleExitValues",
            "RuleBalanceGatheringFlows",
            "EqGeneral",
            "js",
            "costpluscurrents"
        };
        missingKeys = Select[requiredKeys, !KeyExistsQ[Eqs, #]&];
        If[missingKeys =!= {},
            report = <|
                "Valid" -> False,
                "Reason" -> "MissingRequiredFields",
                "MissingKeys" -> missingKeys,
                "BlockChecks" -> <||>,
                "EquationResidual" -> Missing["NotAvailable"],
                "EquationResidualPass" -> False,
                "Tolerance" -> tol
            |>;
            If[verboseQ,
                Print["IsCriticalSolution: missing required fields: ", missingKeys]
            ];
            Return[If[returnReportQ, report, False], Module]
        ];

        solution = CriticalSolutionAssociation[Eqs];

        If[!AssociationQ[solution],
            report = <|
                "Valid" -> False,
                "Reason" -> "MissingAssoCritical",
                "MissingKeys" -> {},
                "BlockChecks" -> <||>,
                "EquationResidual" -> Missing["NotAvailable"],
                "EquationResidualPass" -> False,
                "Tolerance" -> tol
            |>;
            If[verboseQ, Print["IsCriticalSolution: missing association in \"AssoCritical\" or \"Solution\"."]];
            Return[If[returnReportQ, report, False], Module]
        ];

        flowApprox = CriticalFlowApproximation[Eqs, solution];
        costpluscurrents = Lookup[Eqs, "costpluscurrents", Missing["NotAvailable"]];
        eqGeneralCritical =
            Lookup[Eqs, "EqGeneral", True] /.
                If[AssociationQ[costpluscurrents],
                    Expand /@ (costpluscurrents /. flowApprox),
                    {}
                ];

        blocks = <|
            "EqEntryIn" -> (And @@ Lookup[Eqs, "EqEntryIn", True]),
            "EqValueAuxiliaryEdges" -> Lookup[Eqs, "EqValueAuxiliaryEdges", True],
            "AltOptCond" -> Lookup[Eqs, "AltOptCond", True],
            "EqBalanceSplittingFlows" -> Lookup[Eqs, "EqBalanceSplittingFlows", True],
            "AltFlows" -> Lookup[Eqs, "AltFlows", True],
            "AltTransitionFlows" -> Lookup[Eqs, "AltTransitionFlows", True],
            "IneqJs" -> Lookup[Eqs, "IneqJs", True],
            "IneqJts" -> Lookup[Eqs, "IneqJts", True],
            "IneqSwitchingByVertex" -> (And @@ Lookup[Eqs, "IneqSwitchingByVertex", True]),
            "RuleEntryOut" -> RulesConjunction[Lookup[Eqs, "RuleEntryOut", <||>]],
            "RuleExitFlowsIn" -> RulesConjunction[Lookup[Eqs, "RuleExitFlowsIn", <||>]],
            "RuleEntryValues" -> RulesConjunction[Lookup[Eqs, "RuleEntryValues", <||>]],
            "RuleExitValues" -> RulesConjunction[Lookup[Eqs, "RuleExitValues", <||>]],
            "RuleBalanceGatheringFlows" -> RulesConjunction[Lookup[Eqs, "RuleBalanceGatheringFlows", <||>]],
            "EqGeneralCritical" -> eqGeneralCritical
        |>;

        blockResults = Association @ KeyValueMap[
            #1 -> LogicalSatisfiedQ[#2, solution, tol] &,
            blocks
        ];
        concretelyFailed    = Keys @ Select[blockResults, (# === False) &];
        indeterminateBlocks = Keys @ Select[blockResults, (# === Indeterminate) &];
        hasIndeterminate    = indeterminateBlocks =!= {};
        indeterminateResiduals = If[hasIndeterminate,
            KeyTake[Map[Simplify[# /. solution] &, blocks], indeterminateBlocks],
            <||>
        ];

        eqResidualPairs =
            Cases[
                (eqGeneralCritical /. HoldPattern[Equal[lhs_, rhs_]] :> Inactive[Equal][lhs, rhs]) /. solution,
                HoldPattern[Inactive[Equal][lhs_, rhs_]] :> {lhs, rhs},
                Infinity
            ];
        eqResidualDeltas =
            Quiet @ Check[
                N[(#[[1]] - #[[2]])& /@ eqResidualPairs],
                $Failed
            ];
        residual =
            If[eqResidualPairs === {},
                Missing["NotAvailable"],
                If[eqResidualDeltas === $Failed || !VectorQ[eqResidualDeltas, NumericQ],
                    Missing["ComputeError"],
                    Max[Abs[eqResidualDeltas]]
                ]
            ];
        residualPass =
            Which[
                NumericQ[residual],                  residual <= tol,
                residual === Missing["NotAvailable"], TrueQ[Lookup[blockResults, "EqGeneralCritical", False]],
                True,                                Indeterminate
            ];

        overall = Which[
            concretelyFailed =!= {} || residualPass === False, False,
            hasIndeterminate || residualPass === Indeterminate,  Indeterminate,
            True,                                                True
        ];

        If[verboseQ,
            Print["IsCriticalSolution: overall = ", overall];
            If[concretelyFailed =!= {},
                Print["  Concretely failed blocks: ", concretelyFailed]
            ];
            If[indeterminateBlocks =!= {},
                Print["  Indeterminate blocks: ", indeterminateBlocks]
            ];
            Print["IsCriticalSolution: equation residual = ", residual, " (tol=", tol, ")"];
            Print["IsCriticalSolution: equation residual pass = ", residualPass];
        ];

        report = <|
            "Valid"                  -> overall,
            "Reason"                 -> Which[
                                            overall === True,  None,
                                            hasIndeterminate ||
                                            residualPass === Indeterminate, "IndeterminateConstraints",
                                            True,             "ConstraintViolation"],
            "IndeterminateResiduals" -> indeterminateResiduals,
            "ConcretelyFailedBlocks" -> concretelyFailed,
            "MissingKeys"            -> {},
            "BlockChecks"            -> blockResults,
            "EquationResidual"       -> residual,
            "EquationResidualPass"   -> residualPass,
            "Tolerance"              -> tol
        |>;

        (* Backward-compatible boolean: True for valid AND indeterminate, False only for concrete failure *)
        If[returnReportQ, report, !MemberQ[{False}, overall]]
    ];


End[];

EndPackage[];
