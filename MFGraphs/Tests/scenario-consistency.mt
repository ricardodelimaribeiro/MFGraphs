(* Tests for scenario-level consistency checks (switching costs, boundary values). *)

(* Case 12: Core 4x4-type case (Attraction network) - consistent *)
$entries12 = {{1, 100}};
$exits12   = {{4, 0}};

(* Inconsistent Y shortcut - inconsistent *)
$entriesY  = {{1, 100}};
$exitsY    = {{3, 0}, {4, 0}};

Test[
    Module[{s, sys},
        s = getExampleScenario[12, $entries12, $exits12, {}];
        sys = makeSystem[s];
        systemData[sys, "ConsistentCosts"] === True
    ]
    ,
    True
    ,
    TestID -> "Scenario consistency: Case 12 is consistent"
]

Test[
    Module[{s, sys},
        s = getExampleScenario["Inconsistent Y shortcut", $entriesY, $exitsY, Automatic];
        sys = makeSystem[s];
        systemData[sys, "ConsistentCosts"] === True
    ]
    ,
    True
    ,
    TestID -> "Scenario consistency: Inconsistent Y shortcut is enforced consistent"
]

Test[
    Module[{s, sys},
        (* Chain with two exits + manual switching cost that was previously flagged incorrectly *)
        s = getExampleScenario["chain with two exits", {{1, 100}}, {{2, 0}, {3, 20}}, {{1, 2, 3, 10}}];
        sys = makeSystem[s];
        systemData[sys, "ConsistentCosts"] === True
    ]
    ,
    True
    ,
    TestID -> "Scenario consistency: Chain with two exits is consistent (bug fix validation)"
]

Test[
    Module[{s, result},
        s = getExampleScenario["Inconsistent Y shortcut", $entriesY, $exitsY, Automatic];
        result = solveScenario[s];
        ListQ[result] && result =!= {}
    ]
    ,
    True
    ,
    TestID -> "Scenario consistency: Inconsistent Y shortcut solver execution (enforced consistent)"
]

(* Test[
    Module[{s, sys, result},
        (* HRF Scenario 1 (Paper case) *)
        s = getExampleScenario["HRF Scenario 1", {{1, 100}, {9, 100}}, {{8, 0}, {10, 0}}, {}];
        sys = makeSystem[s];
        result = solveScenario[s];
        AssociationQ[result] && result["Rules"] =!= {}
    ]
    ,
    True
    ,
    TestID -> "Scenario consistency: HRF Scenario 1 paper case migrated"
] *)

(* Structural invariants relied on by getKirchhoffLinearSystem and makeSystem
   after the dead defensive guards were removed. If a future scenario builder
   ever fails to populate one of these keys, these tests will catch it before
   downstream Join / lookups blow up. Runs against two scenario shapes so
   builder-specific regressions are more likely to surface. *)

Do[
    Module[{label, s, sys},
        {label, s} = pair;
        sys = makeSystem[s];
        Test[ListQ[systemData[sys, "EqEntryIn"]],             True, TestID -> "System invariant: EqEntryIn is a List ("             <> label <> ")"];
        Test[ListQ[systemData[sys, "BalanceGatheringFlows"]], True, TestID -> "System invariant: BalanceGatheringFlows is a List (" <> label <> ")"];
        Test[ListQ[systemData[sys, "BalanceSplittingFlows"]], True, TestID -> "System invariant: BalanceSplittingFlows is a List (" <> label <> ")"];
        Test[AssociationQ[systemData[sys, "RuleEntryIn"]],    True, TestID -> "System invariant: RuleEntryIn is an Association ("   <> label <> ")"];
        Test[AssociationQ[systemData[sys, "RuleEntryOut"]],   True, TestID -> "System invariant: RuleEntryOut is an Association ("  <> label <> ")"];
        Test[AssociationQ[scenarioData[s, "Topology"]],       True, TestID -> "Scenario invariant: Topology is Association ("       <> label <> ")"]
    ],
    {pair, {
        {"case-12",         getExampleScenario[12, $entries12, $exits12, {}]},
        {"chain-two-exits", getExampleScenario["chain with two exits", {{1, 100}}, {{2, 0}, {3, 20}}, {{1, 2, 3, 10}}]}
    }}
];

(* Boundary-aware symmetry folding (Stage 2 of solver experiments).
   These tests pin down: detection on a known-symmetric topology,
   defensive no-op on a topology with no nontrivial automorphisms,
   that addSymmetryEqualities is a no-op when no perms survive, and
   end-to-end that folding cracks case 21 (which times out without it). *)

Test[
    Length[boundaryPreservingAutomorphisms[
        getExampleScenario[21, {{1, 50}, {2, 50}}, {{10, 0}, {11, 0}, {12, 0}}]]] > 0
    ,
    True
    ,
    TestID -> "Symmetry: case 21 with symmetric boundary detects (1<->2)(3<->4)..."
];

Test[
    boundaryPreservingAutomorphisms[
        getExampleScenario[21, {{1, 30}, {2, 70}}, {{10, 0}, {11, 0}, {12, 0}}]]
    ,
    {}
    ,
    TestID -> "Symmetry: case 21 with asymmetric entries returns no perms"
];

Test[
    boundaryPreservingAutomorphisms[
        getExampleScenario["chain with two exits", {{1, 100}}, {{2, 0}, {3, 20}}, {{1, 2, 3, 10}}]]
    ,
    {}
    ,
    TestID -> "Symmetry: chain with asymmetric exits returns no perms"
];

Module[{s, sys},
    s   = getExampleScenario[21, {{1, 30}, {2, 70}}, {{10, 0}, {11, 0}, {12, 0}}];
    sys = makeSystem[s];
    Test[
        addSymmetryEqualities[sys, s] === sys
        ,
        True
        ,
        TestID -> "Symmetry: addSymmetryEqualities is identity when no perms"
    ]
];

Module[{s, sys, sym, res},
    s   = getExampleScenario[21, {{1, 50}, {2, 50}}, {{10, 0}, {11, 0}, {12, 0}}];
    sys = makeSystem[s];
    sym = addSymmetryEqualities[sys, s];
    res = TimeConstrained[activeSetReduceSystem[sym], 30, $TimedOut];
    Test[
        AssociationQ[res] && isValidSystemSolution[sym, res] === True
        ,
        True
        ,
        TestID -> "Symmetry: case 21 with folding solves and validates within 30s"
    ]
];

(* addOracleEqualities is symbolic-only: tests don't need the optional
   numericOracle subpackage. Verify identity on empty Inactive, identity
   on non-converged oracle, and pin behavior on a synthetic oracle
   association that we know is consistent with the system. *)

Module[{s, sys},
    s   = getExampleScenario[12, $entries12, $exits12, {}];
    sys = makeSystem[s];
    Test[
        addOracleEqualities[sys, <|"Converged" -> True, "Inactive" -> {}|>] === sys
        ,
        True
        ,
        TestID -> "Oracle: addOracleEqualities is identity on empty Inactive"
    ]
];

Module[{s, sys},
    s   = getExampleScenario[12, $entries12, $exits12, {}];
    sys = makeSystem[s];
    Test[
        addOracleEqualities[sys, <|"Converged" -> False, "Inactive" -> {j[1, 2]}|>] === sys
        ,
        True
        ,
        TestID -> "Oracle: addOracleEqualities is identity when Converged is False"
    ]
];

(* Note: addOracleEqualities runs a FindInstance probe on the linear part
   only -- it can detect over-pruning that violates mass balance, but
   complementarity-level over-pruning (oracle's pinned vars conflict with
   Or-disjunctions) shows up downstream as Residual -> False from the
   active-set solver. Callers are expected to check Lookup[result,
   "Residual"] and fall back to the unpruned system if needed. *)

(* Per-disjunct counter conservation. The new BranchesStartedByDisjunct /
   BranchesKeptByDisjunct counters introduced in PR 1 of the
   disjunct-influence study must sum exactly to their flat counterparts;
   BranchesDroppedByDisjunct may be off by the count of non-Or-branch
   drops (the And-block elem===False case in dnfReduceInstrumentedAnd). *)

Module[{sys, rep, sm},
    sys = makeSystem[gridScenario[{3, 3}, {{1, 100}}, {{9, 0}}]];
    rep = dnfReduceDiagnosticReport[sys, "Timeout" -> 30];
    sm  = rep["Summary"];
    Test[
        sm["BranchesStarted"] === Total[Values[sm["BranchesStartedByDisjunct"]]] &&
        sm["BranchesKept"]    === Total[Values[sm["BranchesKeptByDisjunct"]]] &&
        sm["BranchesDropped"] >= Total[Values[sm["BranchesDroppedByDisjunct"]]]
        ,
        True
        ,
        TestID -> "DisjunctCounters: conservation on Grid0303"
    ]
];

Module[{sys, rep, sm},
    sys = makeSystem[getExampleScenario[22, {{1, 100}}, {{6, 0}, {7, 0}}]];
    rep = dnfReduceDiagnosticReport[sys, "Timeout" -> 60];
    sm  = rep["Summary"];
    Test[
        sm["BranchesStarted"] === Total[Values[sm["BranchesStartedByDisjunct"]]] &&
        sm["BranchesKept"]    === Total[Values[sm["BranchesKeptByDisjunct"]]] &&
        sm["BranchesDropped"] >= Total[Values[sm["BranchesDroppedByDisjunct"]]]
        ,
        True
        ,
        TestID -> "DisjunctCounters: conservation on case 22"
    ]
];

(* PR 3 of disjunct-influence study: each Block-* DisjunctOrdering must
   produce the same solution as the Lexicographic baseline. The orderings
   only permute disjuncts, so the rule set (after KeySort) must be identical. *)
sortedRules[res_] := Sort[Lookup[res, "Rules", {}]];

Module[{sys, base},
    sys  = makeSystem[getExampleScenario[12, {{1, 100}}, {{4, 0}}]];
    base = sortedRules[activeSetReduceSystem[sys]];
    Test[
        sortedRules[activeSetReduceSystem[sys, "DisjunctOrdering" -> "Block-Vertex"]],
        base,
        TestID -> "DisjunctOrdering: Block-Vertex matches Lexicographic on case 12"
    ];
    Test[
        sortedRules[activeSetReduceSystem[sys, "DisjunctOrdering" -> "Block-Edge"]],
        base,
        TestID -> "DisjunctOrdering: Block-Edge matches Lexicographic on case 12"
    ];
    Test[
        sortedRules[activeSetReduceSystem[sys, "DisjunctOrdering" -> "Block-SCC"]],
        base,
        TestID -> "DisjunctOrdering: Block-SCC matches Lexicographic on case 12"
    ]
];

Module[{sys, base},
    sys  = makeSystem[getExampleScenario["Grid0303", Automatic, Automatic]];
    base = sortedRules[activeSetReduceSystem[sys]];
    Test[
        sortedRules[activeSetReduceSystem[sys, "DisjunctOrdering" -> "Block-Vertex"]],
        base,
        TestID -> "DisjunctOrdering: Block-Vertex matches Lexicographic on Grid0303"
    ];
    Test[
        sortedRules[activeSetReduceSystem[sys, "DisjunctOrdering" -> "Block-Edge"]],
        base,
        TestID -> "DisjunctOrdering: Block-Edge matches Lexicographic on Grid0303"
    ];
    Test[
        sortedRules[activeSetReduceSystem[sys, "DisjunctOrdering" -> "Block-SCC"]],
        base,
        TestID -> "DisjunctOrdering: Block-SCC matches Lexicographic on Grid0303"
    ]
];
