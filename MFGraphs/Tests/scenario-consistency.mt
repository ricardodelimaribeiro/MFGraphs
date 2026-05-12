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
