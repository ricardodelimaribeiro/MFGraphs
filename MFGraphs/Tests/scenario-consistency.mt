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
        sys = Quiet[makeSystem[s], {mfgSystem::switchingcosts}];
        systemData[sys, "ConsistentCosts"] === False
    ]
    ,
    True
    ,
    TestID -> "Scenario consistency: Inconsistent Y shortcut is detected"
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
    Module[{s, sys, result},
        s = getExampleScenario["Inconsistent Y shortcut", $entriesY, $exitsY, Automatic];
        sys = Quiet[makeSystem[s], {mfgSystem::switchingcosts}];
        result = reduceSystem[sys];
        (* Solver should still return a result, even if inconsistent, but this specific case is feasible *)
        AssociationQ[result] && result["Rules"] =!= {}
    ]
    ,
    True
    ,
    TestID -> "Scenario consistency: Inconsistent Y shortcut solver execution"
]

(* Test[
    Module[{s, sys, result},
        (* HRF Scenario 1 (Paper case) *)
        s = getExampleScenario["HRF Scenario 1", {{1, 100}, {9, 100}}, {{8, 0}, {10, 0}}, {}];
        sys = makeSystem[s];
        result = reduceSystem[sys];
        AssociationQ[result] && result["Rules"] =!= {}
    ]
    ,
    True
    ,
    TestID -> "Scenario consistency: HRF Scenario 1 paper case migrated"
] *)
