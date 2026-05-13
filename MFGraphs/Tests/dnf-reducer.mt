(* dnf-reducer.mt — unit tests for harvestDNFBranch and parseDNFReduceResult *)
Needs["MFGraphs`"];

(* Helper: reconstruct a testable Boolean expression from a parsed solver result *)
toExpr[result_List]        := And @@ (Equal @@@ result);
toExpr[result_Association] := And[
    And @@ (Equal @@@ Lookup[result, "Rules", {}]),
    Lookup[result, "Residual", True]
];
toExpr[result_] := result;


(* ------------------------------------------------------------------ *)
(* harvestDNFBranch                                                    *)
(* ------------------------------------------------------------------ *)

(* Inequality-only residual: Variables[x >= 0] == {}, so the old code
   would call Simplify[x >= 0], get a non-True non-False result, and
   return False — killing a live branch. *)
Test[
    solversTools`Private`harvestDNFBranch[x >= 0, {x}],
    <|"Rules" -> {}, "Residual" -> x >= 0|>,
    TestID -> "harvestDNFBranch: inequality-only branch survives"
]

(* Mixed: one promoted rule, one inequality residual. *)
Test[
    solversTools`Private`harvestDNFBranch[(y == 1) && (x >= 0), {x, y}],
    <|"Rules" -> {y -> 1}, "Residual" -> x >= 0|>,
    TestID -> "harvestDNFBranch: rule plus inequality residual"
]

(* Contradictory equalities must still return False. *)
Test[
    solversTools`Private`harvestDNFBranch[x == 1 && x == 2, {x}],
    False,
    TestID -> "harvestDNFBranch: contradictory equalities return False"
]

(* Fully determined equality: clean rules, Equations -> True. *)
Test[
    solversTools`Private`harvestDNFBranch[x == 3, {x}],
    <|"Rules" -> {x -> 3}, "Residual" -> True|>,
    TestID -> "harvestDNFBranch: fully determined equality returns clean rules"
]


(* ------------------------------------------------------------------ *)
(* parseDNFReduceResult — semantic equivalence to BooleanConvert DNF   *)
(* ------------------------------------------------------------------ *)

(* Mixed inequality + equality disjunction. *)
Test[
    TrueQ @ Resolve[
        ForAll[{x, y},
            Equivalent[
                toExpr[solversTools`Private`parseDNFReduceResult[
                    solversTools`dnfReduce[True, (x >= 0) && ((x == 0) || (y == 1))],
                    {x, y}
                ]],
                BooleanConvert[(x >= 0) && ((x == 0) || (y == 1)), "DNF"]
            ]
        ],
        Reals
    ],
    True,
    TestID -> "parseDNFReduceResult: mixed inequality-equality DNF equivalence"
]

(* Common rule extraction must not collapse unrelated inequality branches. *)
Test[
    TrueQ @ Resolve[
        ForAll[{x, y, z},
            Equivalent[
                toExpr[solversTools`Private`parseDNFReduceResult[
                    solversTools`dnfReduce[
                        True,
                        ((x == 1) && (y >= 0)) || ((x == 1) && (z >= 0))
                    ],
                    {x, y, z}
                ]],
                BooleanConvert[((x == 1) && (y >= 0)) || ((x == 1) && (z >= 0)), "DNF"]
            ]
        ],
        Reals
    ],
    True,
    TestID -> "parseDNFReduceResult: common rule extraction preserves inequality branches"
]

(* Regression: pure equality disjunction must continue to work. TautologyQ
   suffices here (no inequality atoms); the other two tests in this block
   need Resolve[..., Reals] because TautologyQ refuses non-Boolean atoms. *)
Test[
    TautologyQ @ Equivalent[
        toExpr[solversTools`Private`parseDNFReduceResult[
            solversTools`dnfReduce[True, (x == 0) || (x == 1)],
            {x}
        ]],
        BooleanConvert[(x == 0) || (x == 1), "DNF"]
    ],
    True,
    TestID -> "parseDNFReduceResult: pure equality Or preserved (regression)"
]


(* ------------------------------------------------------------------ *)
(* bfsDNFReduce — breadth-first analog of dnfReduceProcedural          *)
(* ------------------------------------------------------------------ *)

(* Pure equality Or: same single-conjunct fan-out as DFS. *)
Test[
    TautologyQ @ Equivalent[
        solversTools`bfsDNFReduce[True, (x == 0) || (x == 1), {x}],
        solversTools`dnfReduce[True, (x == 0) || (x == 1)]
    ],
    True,
    TestID -> "bfsDNFReduce: pure equality Or matches dnfReduce"
]

(* Mixed inequality + Or with shared equality: BFS may emit branches in
   a different order, but the disjunction must be semantically equivalent. *)
Test[
    TrueQ @ Resolve[
        ForAll[{x, y, z},
            Equivalent[
                solversTools`bfsDNFReduce[
                    True,
                    (x >= 0) && ((x == 0) || (y == 1)) && ((y == 1) || (z == 2)),
                    {x, y, z}
                ],
                solversTools`dnfReduce[
                    True,
                    (x >= 0) && ((x == 0) || (y == 1)) && ((y == 1) || (z == 2))
                ]
            ]
        ],
        Reals
    ],
    True,
    TestID -> "bfsDNFReduce: nested Or with shared equality matches DFS semantically"
]

(* parseDNFReduceResult round-trip: BFS output must validate via the same
   harvester as DFS. *)
Test[
    With[{
        bfsResult = solversTools`bfsDNFReduce[
            True,
            ((x == 1) && (y >= 0)) || ((x == 1) && (z >= 0)),
            {x, y, z}
        ]
    },
        TrueQ @ Resolve[
            ForAll[{x, y, z},
                Equivalent[
                    toExpr[solversTools`Private`parseDNFReduceResult[bfsResult, {x, y, z}]],
                    BooleanConvert[((x == 1) && (y >= 0)) || ((x == 1) && (z >= 0)), "DNF"]
                ]
            ],
            Reals
        ]
    ],
    True,
    TestID -> "bfsDNFReduce: parseDNFReduceResult round-trip matches BooleanConvert"
]

(* Infeasible system: BFS must return False, same as DFS. *)
Test[
    solversTools`bfsDNFReduce[True, (x == 0) && (x == 1), {x}],
    False,
    TestID -> "bfsDNFReduce: contradictory equalities return False"
]
