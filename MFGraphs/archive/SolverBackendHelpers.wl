(* Wolfram Language Package *)
(*
   ARCHIVED — solver-phase backend helpers.

   These functions were extracted from MFGraphs/MFGraphs.wl during the
   core-scenario-kernel phase cleanup (2026-04-24). They support the
   Monotone and Critical-Quadratic solver backends (FictitiousPlay, Monotone),
   which are not loaded in the current phase.

   To restore: move these definitions back into MFGraphs.wl Private section
   (or into a dedicated solver-backend submodule) and re-expose any public
   symbols via usage declarations in BeginPackage.
*)

BeginPackage["MFGraphs`"];

Begin["`Private`"];

(* --- Numeric predicate --- *)

NumberVectorQ[j_] := VectorQ[j, NumericQ];

(* --- Feasibility helpers --- *)

IsFeasible[assoc_Association] :=
    Lookup[assoc, "Feasibility", "Infeasible"] === "Feasible";

CheckFlowFeasibility[Null] := "Infeasible";
CheckFlowFeasibility[assoc_Association] :=
    Module[{flowKeys, flowVals},
        flowKeys = Select[Keys[assoc], MatchQ[#, _j] &];
        flowVals = Lookup[assoc, flowKeys];
        If[flowVals === {} || Min[Select[flowVals, NumericQ]] < 0,
            "Infeasible", "Feasible"]
    ];

(* --- Monotone backend: Hessian helpers --- *)

Hess[j_?NumberVectorQ] := DiagonalMatrix[1/# & /@ j];

InverseHessian[j_?NumberVectorQ] := DiagonalMatrix[j];

(* --- Lookup helper --- *)

LookupAssociationValue[assoc_Association, key_, default_:0.] :=
  If[KeyExistsQ[assoc, key], assoc[key], default];

(* --- Monotone backend: state data --- *)

BuildMonotoneStateData[d2e_Association] :=
  Module[{B, KM, jj, halfPairs, signedFlows, rules, substituted, S},
    If[!NameQ["MFGraphs`GetKirchhoffLinearSystem"],
      Return[
        <|
          "B" -> Missing["NotAvailable"],
          "KM" -> Missing["NotAvailable"],
          "FullVariables" -> {},
          "HalfPairs" -> List @@@ Lookup[d2e, "edgeList", {}],
          "SignedEdgeMatrix" -> Missing["NotAvailable"]
        |>,
        Module
      ]
    ];
    {B, KM, jj} = GetKirchhoffLinearSystem[d2e];
    halfPairs = List @@@ Lookup[d2e, "edgeList", {}];
    If[Length[jj] === 0 || Length[halfPairs] === 0,
      S = SparseArray[{}, {Length[halfPairs], Length[jj]}],
      signedFlows = Lookup[d2e, "SignedFlows", <||>];
      rules = Join[
        Lookup[d2e, "RuleBalanceGatheringFlows", <||>],
        Lookup[d2e, "RuleExitFlowsIn", <||>],
        Lookup[d2e, "RuleEntryOut", <||>]
      ];
      substituted = (signedFlows[#] & /@ halfPairs) /. rules;
      S = Last @ CoefficientArrays[substituted, jj]
    ];
    <|
      "B" -> B,
      "KM" -> KM,
      "FullVariables" -> jj,
      "HalfPairs" -> halfPairs,
      "SignedEdgeMatrix" -> S
    |>
  ];

(* --- Monotone backend: value system --- *)

BuildMonotoneValueSystem[d2e_Association] :=
  Module[{statePairs, pairIndex, equalityPairs, boundaryValues, transitions, switching,
      stateCount, transitionsAt},
    statePairs = DeleteDuplicates @ Cases[Lookup[d2e, "us", {}], u[a_, b_] :> {a, b}];
    stateCount = Length[statePairs];
    If[stateCount === 0,
      Return[<|"StateValueAssociation" -> (<||> &)|>, Module]
    ];
    pairIndex = AssociationThread[statePairs, Range[stateCount]];
    equalityPairs = Cases[Lookup[d2e, "Equ", {}], (u @@@ p1_) == (u @@@ p2_) :> {p1, p2}];
    boundaryValues = Cases[Lookup[d2e, "Equ", {}], (u @@@ p_) == val_?NumericQ :> {p, val}];
    transitions = Lookup[d2e, "Transitions", <||>];
    switching = Lookup[d2e, "SwitchingCosts", <||>];
    transitionsAt = GroupBy[Lookup[d2e, "auxTriples", {}], #[[{1, 2}]] &];

    <|
      "StateValueAssociation" -> Function[{edgeCosts},
        Module[{vars, eqs, sol, usLocal},
          vars = u @@@ statePairs;
          eqs = Flatten @ {
            Map[
              Function[{pair},
                u @@ pair == (
                  LookupAssociationValue[transitions, pair, 0] +
                  Plus @@ Map[
                    Function[{triple},
                      u @@ triple[[{2, 3}]] + LookupAssociationValue[switching, triple, 0]
                    ],
                    Lookup[transitionsAt, Key[pair], {}]
                  ]
                )
              ],
              statePairs
            ],
            Map[(u @@@ #[[1]] == u @@@ #[[2]]) &, equalityPairs],
            Map[(u @@@ #[[1]] == #[[2]]) &, boundaryValues]
          };
          usLocal = Lookup[d2e, "us", {}];
          sol = Quiet @ Check[Solve[eqs, usLocal], {}];
          If[MatchQ[sol, {{_Rule ..}}],
            Association[First[sol]],
            Missing["NotAvailable"]
          ]
        ]
      ]
    |>
  ];

(* --- Monotone backend: reduced Kirchhoff coordinates --- *)

BuildReducedKirchhoffCoordinates[d2e_Association, basePoint_:Automatic] :=
  Module[{stateData, KM, jj, S, nullBasis, invisibleBasis, visibleBasis,
      removeInvisible, basisMatrix, dim, baseVec, edgeOffset, edgeBasis,
      invisibleDim, signedEdgeRows, tol = 10^-10},
    stateData = BuildMonotoneStateData[d2e];
    KM = stateData["KM"];
    jj = stateData["FullVariables"];
    S = stateData["SignedEdgeMatrix"];
    nullBasis = N @ Orthogonalize[NullSpace[Normal[KM]]];
    invisibleBasis = N @ Orthogonalize[NullSpace[Join[Normal[KM], Normal[S]]]];
    removeInvisible[v_] :=
      Fold[#1 - (#1.#2) #2 &, v, invisibleBasis];
    visibleBasis = If[nullBasis === {},
      {},
      If[invisibleBasis === {},
        nullBasis,
        Orthogonalize @ Select[removeInvisible /@ nullBasis, Norm[#] > tol &]
      ]
    ];
    visibleBasis = Select[visibleBasis, Norm[#] > tol &];
    basisMatrix = If[visibleBasis === {},
      ConstantArray[0., {Length[jj], 0}],
      N @ Transpose[visibleBasis]
    ];
    dim = Length[visibleBasis];
    baseVec = Which[
      basePoint === Automatic && Length[jj] === 0, {},
      basePoint === Automatic, Quiet @ Check[
        N @ LeastSquares[N @ Normal[KM], N @ stateData["B"]],
        Missing["NotAvailable"]
      ],
      True, N @ basePoint
    ];
    signedEdgeRows = If[ArrayQ[S], First @ Dimensions[S], 0];
    edgeOffset = If[signedEdgeRows === 0,
      {},
      If[ListQ[baseVec] && VectorQ[baseVec, NumericQ],
        N @ (S . baseVec),
        Missing["NotAvailable"]
      ]
    ];
    edgeBasis = If[signedEdgeRows === 0,
      ConstantArray[0., {0, dim}],
      If[MatrixQ[basisMatrix], N @ (S . basisMatrix), Missing["NotAvailable"]]
    ];
    invisibleDim = Length[invisibleBasis];
    <|
      "BasePoint" -> baseVec,
      "BasisMatrix" -> basisMatrix,
      "StateDimension" -> dim,
      "FullDimension" -> Length[jj],
      "CostInvisibleDimension" -> invisibleDim,
      "SignedEdgeOffset" -> edgeOffset,
      "SignedEdgeBasis" -> edgeBasis,
      "FullVariables" -> jj
    |>
  ];

(* --- Critical-quadratic backend --- *)

BuildMonotonePairCostAssociation[halfPairs_List, edgeList_List, q_] :=
  Association @ Flatten @ MapThread[
    Function[{pair, edge, flow},
      With[{cost = If[PossibleZeroQ[flow], 0., N @ Abs[flow]]},
        {pair -> cost, Reverse[pair] -> cost}
      ]
    ],
    {halfPairs, edgeList, q}
  ];

BuildCriticalQuadraticEdgeModel[d2e_Association, tol_:10^-8] :=
  Module[{edges},
    edges = Lookup[d2e, "edgeList", {}];
    AssociationThread[edges, ConstantArray[1., Length[edges]]]
  ];

BuildCriticalQuadraticObjective[d2e_Association] :=
  Module[{stateData, B, KM, jj, S, q, edgeSlopes, switching, exitRules,
      outExitPairs, exitCostByPair, linearTerm, quadraticTerm, edgeSlopeVector},
    stateData = BuildMonotoneStateData[d2e];
    B = stateData["B"];
    KM = stateData["KM"];
    jj = stateData["FullVariables"];
    If[Length[jj] === 0,
      Return[$Failed, Module]
    ];
    edgeSlopes = BuildCriticalQuadraticEdgeModel[d2e];
    If[edgeSlopes === $Failed,
      Return[$Failed, Module]
    ];
    S = stateData["SignedEdgeMatrix"];
    q = If[Length[stateData["HalfPairs"]] === 0, {}, S . jj];
    edgeSlopeVector = Lookup[edgeSlopes, Lookup[d2e, "edgeList", {}], 1.];
    quadraticTerm =
      If[q === {} || edgeSlopeVector === {},
        0,
        1/2 q . DiagonalMatrix[edgeSlopeVector] . q
      ];
    switching = Lookup[d2e, "SwitchingCosts", <||>];
    exitRules = Lookup[d2e, "RuleExitValues", <||>];
    outExitPairs = List @@@ Lookup[d2e, "exitEdges", {}];
    exitCostByPair =
      Association @ Table[
        pair -> LookupAssociationValue[exitRules, u @@ pair, 0],
        {pair, Join[outExitPairs, Reverse /@ outExitPairs]}
      ];
    linearTerm =
      Total[
        (
          LookupAssociationValue[switching, #, 0] +
          LookupAssociationValue[exitCostByPair, #[[{2, 3}]], 0]
        ) (j @@ #) & /@ Lookup[d2e, "auxTriples", {}]
      ] +
      10^-7 Total[jj];
    <|
      "B" -> B,
      "KM" -> KM,
      "Variables" -> jj,
      "Objective" -> quadraticTerm + linearTerm,
      "StateData" -> stateData
    |>
  ];

UseQuadraticCriticalBackendQ[d2e_Association] :=
  Module[{edges, graph, degrees, maxDegree, nonzeroSwitching},
    edges = Lookup[d2e, "edgeList", {}];
    graph = Graph[edges];
    degrees = VertexDegree[graph];
    maxDegree = If[degrees === {}, 0, Max[degrees]];
    nonzeroSwitching =
      AnyTrue[
        Values @ Lookup[d2e, "SwitchingCosts", <||>],
        NumericQ[#] && !PossibleZeroQ[#] &
      ];
    TrueQ[AcyclicGraphQ[graph]] && (maxDegree <= 2 || !nonzeroSwitching)
  ];

(* --- Monotone variable field evaluation --- *)

MonotoneVariableFieldValue[var_, values_Association, switching_Association] :=
  Replace[
    var,
    {
      j[r_, i_, w_] :> N @ (
        LookupAssociationValue[switching, {r, i, w}, 0] +
        LookupAssociationValue[values, u[i, w], 0.] +
        LookupAssociationValue[values["PairCosts"], {i, w}, 0.]
      ),
      j[a_, b_] :> N @ (
        LookupAssociationValue[values, u[a, b], 0.] +
        LookupAssociationValue[values["PairCosts"], {a, b}, 0.]
      )
    }
  ];

(* --- Flow vector encode/decode --- *)

EncodeFlowAssociation[numericState_Association, assoc_Association] :=
    Module[{vars, vals},
        vars = Lookup[numericState, "FlowVariables", {}];
        If[vars === {},
            Return[Developer`ToPackedArray[{}], Module]
        ];
        vals = vars /. assoc;
        vals = vals /. {u[___] :> 0., j[___] :> 0., z[___] :> 0.};
        vals = Quiet @ Check[N @ vals, vals];
        If[VectorQ[vals, NumericQ],
            Developer`ToPackedArray @ N @ vals,
            Missing["NotAvailable"]
        ]
    ];

EncodeFlowAssociation[_, _] := Missing["NotAvailable"];

DecodeFlowVector[numericState_Association, flowVec_] :=
    Module[{vars, vals},
        vars = Lookup[numericState, "FlowVariables", {}];
        vals = Quiet @ Check[Developer`ToPackedArray @ N @ flowVec, $Failed];
        If[vals === $Failed || !VectorQ[vals, NumericQ] || Length[vals] =!= Length[vars],
            <||>,
            AssociationThread[vars, vals]
        ]
    ];

DecodeFlowVector[_, _] := <||>;

(* --- Fast Kirchhoff residual / signed-edge flow computation --- *)

ComputeSignedEdgeFlowsFast[numericState_Association, flowVec_] :=
    Module[{indices, signedEdgeMatrix, kirchhoffVec},
        indices = Lookup[numericState, "KirchhoffIndicesInFlowVector", {}];
        signedEdgeMatrix = Lookup[numericState, "SignedEdgeMatrix", SparseArray[{}, {0, 0}]];
        If[Dimensions[signedEdgeMatrix][[1]] === 0,
            Return[Developer`ToPackedArray[{}], Module]
        ];
        If[indices === {} || !VectorQ[indices, IntegerQ],
            Return[Missing["NotAvailable"], Module]
        ];
        If[Length[flowVec] < Max[indices],
            Return[Missing["NotAvailable"], Module]
        ];
        kirchhoffVec = flowVec[[indices]];
        If[!VectorQ[kirchhoffVec, NumericQ],
            Return[Missing["NotAvailable"], Module]
        ];
        Developer`ToPackedArray @ N @ (signedEdgeMatrix . kirchhoffVec)
    ];

ComputeKirchhoffResidualFast[numericState_Association, flowVec_] :=
    Module[{indices, KM, B, kirchhoffVec, residual},
        indices = Lookup[numericState, "KirchhoffIndicesInFlowVector", {}];
        KM = Lookup[numericState, "KirchhoffKM", SparseArray[{}, {0, 0}]];
        B = Lookup[numericState, "KirchhoffB", Developer`ToPackedArray[{}]];
        If[indices === {} || !VectorQ[indices, IntegerQ],
            Return[Missing["NotAvailable"], Module]
        ];
        If[Length[flowVec] < Max[indices],
            Return[Missing["NotAvailable"], Module]
        ];
        kirchhoffVec = flowVec[[indices]];
        If[!VectorQ[kirchhoffVec, NumericQ],
            Return[Missing["NotAvailable"], Module]
        ];
        residual = Quiet @ Check[N @ Norm[KM . kirchhoffVec - B, Infinity], Missing["NotAvailable"]];
        If[NumericQ[residual], residual, Missing["NotAvailable"]]
    ];

(* --- Flow selection --- *)

SelectFlowAssociation[assoc_Association] :=
    Association @ KeySelect[assoc, MatchQ[#, _j] &];

(* --- Boundary mass data --- *)

BuildBoundaryMassData[Eqs_Association, flowAssoc_Association] :=
    Module[{missing, entryPairs, exitPairs, entryFlows, exitFlows, inTotal, outTotal, residual},
        missing = Missing["NotAvailable"];
        entryPairs = List @@@ Lookup[Eqs, "entryEdges", {}];
        exitPairs = List @@@ Lookup[Eqs, "exitEdges", {}];
        entryFlows =
            If[entryPairs === {},
                {},
                (j @@@ entryPairs) /. flowAssoc
            ];
        exitFlows =
            If[exitPairs === {},
                {},
                (j @@@ exitPairs) /. flowAssoc
            ];
        entryFlows = entryFlows /. {j[___] :> 0., u[___] :> 0.};
        exitFlows = exitFlows /. {j[___] :> 0., u[___] :> 0.};
        inTotal =
            If[
                ListQ[entryFlows] && VectorQ[entryFlows, NumericQ],
                N @ Total[entryFlows],
                missing
            ];
        outTotal =
            If[
                ListQ[exitFlows] && VectorQ[exitFlows, NumericQ],
                N @ Total[exitFlows],
                missing
            ];
        residual =
            If[
                NumericQ[inTotal] && NumericQ[outTotal],
                N @ Abs[inTotal - outTotal],
                missing
            ];
        <|
            "BoundaryMassIn" -> inTotal,
            "BoundaryMassOut" -> outTotal,
            "BoundaryMassResidual" -> residual
        |>
    ];

BuildBoundaryMassData[_, _] :=
    <|
        "BoundaryMassIn" -> Missing["NotAvailable"],
        "BoundaryMassOut" -> Missing["NotAvailable"],
        "BoundaryMassResidual" -> Missing["NotAvailable"]
    |>;

(* --- Utility reduction residual data --- *)

BuildUtilityReductionResidualData[Eqs_Association, solution_Association] :=
    Module[{missing, utilityReduction, classes, classResiduals, overallResidual,
      classCount, reducedVarCount},
        missing = Missing["NotAvailable"];
        utilityReduction = Lookup[Eqs, "UtilityReduction", Missing["NotAvailable"]];
        If[!AssociationQ[utilityReduction],
            Return[
                <|
                    "UtilityClassResidual" -> missing,
                    "UtilityClassCount" -> missing,
                    "UtilityReducedVariableCount" -> missing
                |>,
                Module
            ]
        ];
        classes = Lookup[utilityReduction, "UtilityClasses", {}];
        classCount =
            If[ListQ[classes], Length[classes], missing];
        reducedVarCount = Lookup[utilityReduction, "ReducedVariableCount", missing];
        classResiduals = Cases[
            Map[
                Function[class,
                    Module[{vals},
                        vals = Quiet @ Check[N @ (class /. solution), $Failed];
                        If[ListQ[vals] && VectorQ[vals, NumericQ] && Length[vals] > 1,
                            N @ Max @ Abs[vals - First[vals]],
                            If[ListQ[vals] && VectorQ[vals, NumericQ], 0., Missing["NotAvailable"]]
                        ]
                    ]
                ],
                classes
            ],
            _?NumericQ
        ];
        overallResidual =
            If[classResiduals === {}, missing, N @ Max[classResiduals]];
        <|
            "UtilityClassResidual" -> overallResidual,
            "UtilityClassCount" -> classCount,
            "UtilityReducedVariableCount" -> reducedVarCount
        |>
    ];

BuildUtilityReductionResidualData[_, _] :=
    <|
        "UtilityClassResidual" -> Missing["NotAvailable"],
        "UtilityClassCount" -> Missing["NotAvailable"],
        "UtilityReducedVariableCount" -> Missing["NotAvailable"]
    |>;

(* --- Solver comparison data builder --- *)

BuildSolverComparisonData[Eqs_Association, solution_] :=
    Module[{missing, flowAssoc, edgeList, edgePairs, signedFlowRules, signedExprs,
         signedVals, B, KM, jj, jjVals, residual, numericState, encodedFlowVec,
         fastSignedVals, fastResidual, usedFastSignedQ, usedFastResidualQ,
         boundaryMassData, utilityReductionResidualData},
        missing = Missing["NotAvailable"];
        edgeList = Lookup[Eqs, "edgeList", {}];
        If[!AssociationQ[solution],
            Return[
                <|
                    "ComparableEdges" -> edgeList,
                    "FlowAssociation" -> missing,
                    "SignedEdgeFlows" -> missing,
                    "ComparableFlowVector" -> missing,
                    "KirchhoffResidual" -> missing,
                    "BoundaryMassIn" -> missing,
                    "BoundaryMassOut" -> missing,
                    "BoundaryMassResidual" -> missing,
                    "UtilityClassResidual" -> missing,
                    "UtilityClassCount" -> missing,
                    "UtilityReducedVariableCount" -> missing
                |>,
                Module
            ]
        ];
        flowAssoc = SelectFlowAssociation[solution];
        boundaryMassData = BuildBoundaryMassData[Eqs, flowAssoc];
        utilityReductionResidualData = BuildUtilityReductionResidualData[Eqs, solution];
        numericState = Lookup[Eqs, "NumericState", Missing["NotAvailable"]];
        usedFastSignedQ = False;
        usedFastResidualQ = False;
        signedVals = missing;
        residual = missing;

        If[AssociationQ[numericState],
            encodedFlowVec = Quiet @ Check[
                EncodeFlowAssociation[numericState, flowAssoc],
                missing
            ];
            If[ListQ[encodedFlowVec] && VectorQ[encodedFlowVec, NumericQ],
                fastSignedVals = Quiet @ Check[
                    ComputeSignedEdgeFlowsFast[numericState, encodedFlowVec],
                    missing
                ];
                If[
                    ListQ[fastSignedVals] &&
                    VectorQ[fastSignedVals, NumericQ] &&
                    Length[fastSignedVals] === Length[edgeList],
                    signedVals = fastSignedVals;
                    usedFastSignedQ = True
                ];
                fastResidual = Quiet @ Check[
                    ComputeKirchhoffResidualFast[numericState, encodedFlowVec],
                    missing
                ];
                If[NumericQ[fastResidual],
                    residual = fastResidual;
                    usedFastResidualQ = True
                ]
            ]
        ];

        If[!usedFastSignedQ,
            edgePairs = List @@@ edgeList;
            signedFlowRules = Lookup[Eqs, "SignedFlows", <||>];
            signedExprs =
                (If[KeyExistsQ[signedFlowRules, #], signedFlowRules[#], missing] & /@ edgePairs) /.
                    Join[
                        Lookup[Eqs, "RuleBalanceGatheringFlows", <||>],
                        Lookup[Eqs, "RuleExitFlowsIn", <||>],
                        Lookup[Eqs, "RuleEntryOut", <||>]
                    ];
            signedVals = Quiet @ Check[signedExprs /. flowAssoc, missing];
            If[!ListQ[signedVals] || !VectorQ[signedVals, NumericQ],
                signedVals = missing
            ]
        ];

        If[!usedFastResidualQ && NameQ["MFGraphs`GetKirchhoffLinearSystem"],
            {B, KM, jj} = GetKirchhoffLinearSystem[Eqs];
            jjVals = Lookup[flowAssoc, jj, missing];
            residual =
                If[ListQ[jjVals] && VectorQ[jjVals, NumericQ],
                    N @ Norm[KM . jjVals - B, Infinity],
                    missing
                ]
        ];
        Join[
            <|
                "ComparableEdges" -> edgeList,
                "FlowAssociation" -> flowAssoc,
                "SignedEdgeFlows" ->
                    If[signedVals === missing, missing, AssociationThread[edgeList, signedVals]],
                "ComparableFlowVector" -> signedVals,
                "KirchhoffResidual" -> residual
            |>,
            boundaryMassData,
            utilityReductionResidualData
        ]
    ];

(* --- Solver result constructor --- *)

MakeSolverResult[
    solver_String,
    resultKind_String,
    feasibility_,
    message_,
    solution_,
    extraAssoc_: <||>
] :=
    Join[
        <|
            "Solver" -> solver,
            "ResultKind" -> resultKind,
            "Feasibility" -> feasibility,
            "Message" -> message,
            "Solution" -> solution
        |>,
        <|"Convergence" -> Missing["NotApplicable"]|>,
        extraAssoc
    ];

(* --- Compiled input predicate --- *)

SolveMFGCompiledInputQ[input_Association] :=
    KeyExistsQ[input, "EqGeneral"] &&
    KeyExistsQ[input, "js"] &&
    KeyExistsQ[input, "jts"];

(* --- ReindexToIntegers — archived with review findings (2026-04-25) ---

   DECISION: makeScenario hard-rejects non-integer vertex labels (returns Failure).
   This function was implemented and then archived. Re-enable if soft-reindex is
   ever preferred over hard-reject at the API boundary.

   Known limitations of this implementation:
   - Association-form switching costs pass straight through (no key remapping). If the
     caller supplies <|{i,k,j}->cost|> with non-integer keys, those keys are not remapped.
   - Reindexing runs before makeScenario validation; a bad Graph argument could panic inside
     VertexList/AdjacencyMatrix before the Failure path is reached.
   - Does not handle the "Graph" key path: scenarioFromGraph now passes the Graph object
     directly and relies on NormalizeScenarioModel to derive Vertices/Adjacency.
*)
ReindexToIntegers[vl_, entries_, exits_, sc_] :=
    If[AllTrue[vl, IntegerQ],
        {vl, entries, exits, sc},
        Module[{idx},
            idx = AssociationThread[vl, Range[Length[vl]]];
            {
                Range[Length[vl]],
                {idx[First[#]], Last[#]} & /@ entries,
                {idx[First[#]], Last[#]} & /@ exits,
                If[ListQ[sc] && sc =!= {},
                    {idx[#[[1]]], idx[#[[2]]], idx[#[[3]]], #[[4]]} & /@ sc,
                    sc
                ]
            }
        ]
    ];

End[];

EndPackage[];
