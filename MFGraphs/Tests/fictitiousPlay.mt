(* Tests for the K-fold consensus oracle in numericOracle.wl.
   OPT-IN: invokes Needs["numericOracle`"] which loads the float-isolated
   subpackage. Not run by RunTests fast; invoke explicitly via
   RunSingleTest.wls or a future "oracle" tag. *)

Needs["numericOracle`"];

(* fp-converges-grid-3x3: K-fold consensus succeeds and the pruned
   system solves correctly. *)
Module[{s, sys, oracle, pruned, res},
    s      = gridScenario[{3, 3}, {{1, 100}}, {{9, 0}}];
    sys    = makeSystem[s];
    oracle = numericFictitiousPlayClassify[sys,
        "MaxIterations" -> 5, "Timeout" -> 30];
    pruned = addOracleEqualities[sys, oracle];
    res    = TimeConstrained[activeSetReduceSystem[pruned], 30, $TimedOut];
    Test[
        oracle["Converged"] === True &&
        AssociationQ[res] &&
        Lookup[res, "Residual", True] =!= False
        ,
        True
        ,
        TestID -> "FP: grid-3x3 converges and pruned solve is feasible"
    ]
];

(* fp-result-shape-matches-numericOracleClassify: addOracleEqualities
   consumes the FP result without a code path change. *)
Module[{s, sys, fpResult},
    s        = gridScenario[{3, 3}, {{1, 100}}, {{9, 0}}];
    sys      = makeSystem[s];
    fpResult = numericFictitiousPlayClassify[sys,
        "MaxIterations" -> 3, "Timeout" -> 10];
    Test[
        AssociationQ[fpResult] &&
        SubsetQ[Keys[fpResult],
            {"Inactive", "Active", "Ambiguous", "Converged", "Method"}]
        ,
        True
        ,
        TestID -> "FP: result Association has the same keys as numericOracleClassify"
    ]
];

(* fp-cracks-hrf: HRF Scenario 1 cracks via the FP consensus oracle.
   Was INFEASIBLE under the one-shot oracle. *)
Module[{s, sys, oracle, pruned, res},
    s      = getExampleScenario["HRF Scenario 1", {{1, 100}}, {{8, 0}, {10, 0}}];
    sys    = makeSystem[s];
    oracle = numericFictitiousPlayClassify[sys,
        "MaxIterations" -> 5, "Timeout" -> 30];
    pruned = addOracleEqualities[sys, oracle];
    res    = TimeConstrained[activeSetReduceSystem[pruned], 60, $TimedOut];
    Test[
        AssociationQ[res] &&
        res =!= $TimedOut &&
        Lookup[res, "Residual", True] =!= False &&
        isValidSystemSolution[pruned, res] === True
        ,
        True
        ,
        TestID -> "FP: HRF Scenario 1 cracks via consensus + pruning + symbolic solve"
    ]
];
