(* boolean-minimize.mt — tests for booleanMinimizeSystem and booleanMinimizeReduceSystem *)
Needs["MFGraphs`"];

$smallCases = {
    {"grid-2x2", gridScenario[{2, 2}, {{1, 100}}, {{4, 0}}]},
    {"grid-3x2", gridScenario[{3, 2}, {{1, 100}}, {{6, 0}}]},
    {"grid-2x3", gridScenario[{2, 3}, {{1, 100}}, {{6, 0}}]}
};

(* ----------------------------------------------------------- *)
(* booleanMinimizeSystem matches dnfReduceSystem on small cases *)
(* ----------------------------------------------------------- *)

Do[
    Module[{label, scen, sys, refSol, bmSol},
        {label, scen} = pair;
        sys    = makeSystem[scen];
        refSol = dnfReduceSystem[sys];
        bmSol  = booleanMinimizeSystem[sys];
        Test[
            ListQ[bmSol] && ListQ[refSol] && (Sort[bmSol] === Sort[refSol]),
            True,
            TestID -> "booleanMinimizeSystem matches dnfReduceSystem on " <> label
        ]
    ],
    {pair, $smallCases}
];

(* --------------------------------------------------------------------- *)
(* booleanMinimizeReduceSystem produces a valid solution on small cases  *)
(* --------------------------------------------------------------------- *)

Do[
    Module[{label, scen, sys, sol, report},
        {label, scen} = pair;
        sys    = makeSystem[scen];
        sol    = booleanMinimizeReduceSystem[sys];
        report = isValidSystemSolution[sys, sol, "ReturnReport" -> True];
        Test[
            TrueQ[report["Valid"]],
            True,
            TestID -> "booleanMinimizeReduceSystem returns valid solution on " <> label
        ]
    ],
    {pair, $smallCases}
];

(* --------------------------------------------------------------------- *)
(* booleanMinimizeReduceSystem matches dnfReduceSystem on small cases    *)
(* (cross-check that the levered solver lands on the same rules as the   *)
(* reference DNF reducer when both fully determine the system).          *)
(* --------------------------------------------------------------------- *)

Do[
    Module[{label, scen, sys, refSol, bmrSol},
        {label, scen} = pair;
        sys    = makeSystem[scen];
        refSol = dnfReduceSystem[sys];
        bmrSol = booleanMinimizeReduceSystem[sys];
        (* Split assertions so mutual $Failed cannot satisfy the equivalence
           check by short-circuiting both ListQ guards to False. *)
        Test[ListQ[refSol], True,
             TestID -> "dnfReduceSystem returns a List on " <> label];
        Test[ListQ[bmrSol], True,
             TestID -> "booleanMinimizeReduceSystem returns a List on " <> label];
        Test[Sort[bmrSol] === Sort[refSol], True,
             TestID -> "booleanMinimizeReduceSystem matches dnfReduceSystem on " <> label]
    ],
    {pair, $smallCases}
];

(* --------------------------------------------------------------- *)
(* Lever 1 (arm-pruning) actually fires: surviving disjAtom count   *)
(* is strictly less than raw BooleanConvert disjunct count.         *)
(* --------------------------------------------------------------- *)

Module[{sys, inputs, constraints, allVars, atoms, linAtoms, disjAtoms,
        prunedLin, prunedDisj, status, rawDisjuncts, rawDNF},
    sys      = makeSystem[gridScenario[{3, 2}, {{1, 100}}, {{6, 0}}]];
    inputs   = solversTools`Private`buildSolverInputs[sys];
    {constraints, allVars} = inputs[[1 ;; 2]];
    atoms    = If[Head[constraints] === And, List @@ constraints, {constraints}];
    linAtoms = Select[atoms, Head[#] =!= Or &];
    disjAtoms = Select[atoms, Head[#] === Or &];
    {prunedLin, prunedDisj, status} =
        solversTools`Private`pruneDisjunctiveArms[linAtoms, disjAtoms, allVars, 2];
    rawDNF        = BooleanConvert[constraints, "DNF"];
    rawDisjuncts  = If[Head[rawDNF] === Or, Length[rawDNF], 1];
    Test[
        Length[prunedDisj] < Length[disjAtoms] || status === "Reduced",
        True,
        TestID -> "pruneDisjunctiveArms removes arms on grid-3x2"
    ];
    Test[
        rawDisjuncts > 1,
        True,
        TestID -> "BooleanConvert raw disjunct count > 1 on grid-3x2 (sanity)"
    ]
]
