(* Tests for reduceSystem *)

Test[
    NameQ["solversTools`reduceSystem"],
    True,
    TestID -> "reduceSystem: public symbol exists"
]

Test[
    Module[{data, s, sys, result},
        data = <|
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> {}
        |>;
        s = makeScenario[<|"Model" -> data|>];
        sys = makeSystem[s];
        result = reduceSystem[sys];
        result =!= False
    ],
    True,
    TestID -> "reduceSystem: chain 1-exit yields non-False solution"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = reduceSystem[sys];
        result =!= False
    ],
    True,
    TestID -> "reduceSystem: chain 2-exits yields non-False solution"
]

Test[
    Module[{s, sys, entryVals, exitVals},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        entryVals = Values @ Normal @ systemData[sys, "RuleEntryIn"];
        exitVals = Values @ Normal @ systemData[sys, "RuleExitValues"];
        FreeQ[Join[entryVals, exitVals], _Real]
    ],
    True,
    TestID -> "reduceSystem: boundary rules are exactified (no Real coefficients)"
]

Test[
    Module[{s, sys, exitRules},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        exitRules = systemData[sys, "RuleExitValues"];
        KeyExistsQ[exitRules, u["auxExit2", 2]] &&
        KeyExistsQ[exitRules, u["auxExit3", 3]] &&
        !KeyExistsQ[exitRules, u[2, "auxExit2"]] &&
        !KeyExistsQ[exitRules, u[3, "auxExit3"]]
    ],
    True,
    TestID -> "reduceSystem: RuleExitValues use u[auxExit, vertex] orientation"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 10.0}, {3, 0.0}}];
        sys = makeSystem[s];
        result = reduceSystem[sys];
        result =!= False
    ],
    True,
    TestID -> "reduceSystem: chain with costs {2,10},{3,0} remains solvable"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = reduceSystem[sys];
        result =!= False
    ],
    True,
    TestID -> "reduceSystem: chain with costs {2,0},{3,10} remains solvable via inequality+complementarity"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 10.5}}, {{3, 0.25}}];
        sys = makeSystem[s];
        result = reduceSystem[sys];
        result =!= False
    ],
    True,
    TestID -> "reduceSystem: non-integer decimal boundaries solve after exactification"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2];
        sys = makeSystem[s];
        result = Quiet[reduceSystem[sys], reduceSystem::noncritical];
        FailureQ[result] && result["Tag"] === "reduceSystem"
    ],
    True,
    TestID -> "reduceSystem: non-critical congestion systems fail"
]

(* --- dnfReduceSystem tests --- *)

Test[
    NameQ["solversTools`dnfReduceSystem"],
    True,
    TestID -> "dnfReduceSystem: public symbol exists"
]

Test[
    solversTools`Private`branchStateReduceResult[x == 1 && x == 2, {x}],
    <|"Rules" -> {}, "Equations" -> False|>,
    TestID -> "branchStateReduceResult: rejects direct conflicting equalities"
]

Test[
    solversTools`Private`branchStateReduceResult[(x == 1 || x == 2) && x == 3, {x}],
    <|"Rules" -> {}, "Equations" -> False|>,
    TestID -> "branchStateReduceResult: rejects disjuncts conflicting with later equality"
]

Test[
    solversTools`Private`branchStateReduceResult[(x == 1 || x == 2) && x == 2, {x}],
    {x -> 2},
    TestID -> "branchStateReduceResult: keeps compatible disjunct matching later equality"
]

Test[
    Module[{result},
        result = solversTools`Private`branchStateReduceResult[x + y == 1, {x, y}];
        AssociationQ[result] &&
        Lookup[result, "Rules", Missing["Rules"]] === {} &&
        !FreeQ[Lookup[result, "Equations", True], x + y == 1]
    ],
    True,
    TestID -> "branchStateReduceResult: symbolic RHS remains residual"
]

Test[
    solversTools`Private`solutionResultKind[
        solversTools`Private`branchStateReduceResult[x + y == 1, {x, y}]
    ],
    "Residual",
    TestID -> "solutionResultKind: symbolic RHS branch-state result is residual"
]

Test[
    solversTools`Private`branchStateReduceResult[x == y && y == 3, {x, y}],
    {x -> 3, y -> 3},
    TestID -> "branchStateReduceResult: chained equalities become ground rules"
]

Test[
    solversTools`Private`mergeRules[{x -> 1}, {x -> 2}],
    {x -> 2},
    TestID -> "mergeRules: duplicate lhs keeps latest rule"
]

Test[
    solversTools`Private`solutionResultKind[{j[1, 2] -> 10}],
    "Rules",
    TestID -> "solutionResultKind: rule list is Rules"
]

Test[
    solversTools`Private`solutionResultKind[<|"Rules" -> {j[1, 2] -> 10}, "Equations" -> True|>],
    "Rules",
    TestID -> "solutionResultKind: true residual association is Rules"
]

Test[
    solversTools`Private`solutionResultKind[<|"Rules" -> {}, "Equations" -> False|>],
    "NoSolution",
    TestID -> "solutionResultKind: false residual association is NoSolution"
]

Test[
    solversTools`Private`solutionResultKind[<|"Rules" -> {}, "Equations" -> (j[1, 2] == 0 || j[1, 2] == 1)|>],
    "Branched",
    TestID -> "solutionResultKind: residual with Or is Branched"
]

Test[
    solversTools`Private`solutionResultKind[<|"Rules" -> {}, "Equations" -> (j[1, 2] >= 0 && u[1, 2] <= 1)|>],
    "Parametric",
    TestID -> "solutionResultKind: tracked residual without Or is Parametric"
]

Test[
    Module[{diag},
        diag = solversTools`Private`dnfResidualDiagnostics[
            <|"Rules" -> {}, "Equations" -> ((j[1, 2] == 0 && (u[1, 2] == 0 || u[1, 2] == 1)) || j[1, 2] == 1)|>
        ];
        diag["TopLevelBranchCount"] === 2 &&
        diag["NestedOrQ"] === True &&
        diag["TrackedVariablesQ"] === True
    ],
    True,
    TestID -> "dnfResidualDiagnostics: reports branches nested Or and tracked variables"
]

Test[
    Module[{expr, ordered},
        expr = (u[1, 2] >= 0) && ((j[1, 2] == 0) || (u[1, 2] == 0)) && (j[1, 2] == 3);
        ordered = solversTools`Private`dnfOrderConjuncts[expr, "original"];
        ordered === List @@ expr
    ],
    True,
    TestID -> "dnf ordering: original preserves current top-level order"
]

Test[
    Module[{expr, original, ordered},
        expr = (u[1, 2] >= 0) && ((j[1, 2] == 0) || (u[1, 2] == 0)) && (j[1, 2] == 3);
        original = ToString[#, InputForm] & /@ solversTools`Private`dnfTopLevelConjuncts[expr];
        ordered = ToString[#, InputForm] & /@ solversTools`Private`dnfOrderConjuncts[expr, "nonor-stable"];
        Sort[original] === Sort[ordered] && First[ordered] === ToString[j[1, 2] == 3, InputForm]
    ],
    True,
    TestID -> "dnf ordering: nonor-stable preserves conjunct multiset"
]

Test[
    solversTools`Private`dnfVarVertices /@ {j[1, 2], j[1, 2, 3], u[2, 3]},
    {{1, 2}, {1, 2, 3}, {2, 3}},
    TestID -> "dnf ordering: tracked variable vertex extraction"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        (ListQ[result] && MatchQ[result, {__Rule}]) ||
        (AssociationQ[result] && KeyExistsQ[result, "Rules"])
    ],
    True,
    TestID -> "dnfReduceSystem: chain 2-exits returns rules or rules+equations"
]

Test[
    Module[{s, sys, expected, report},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        sys = makeSystem[s];
        expected = dnfReduceSystem[sys];
        report = solversTools`Private`dnfReduceDiagnosticReport[sys, "Order" -> "original", "Timeout" -> 30];
        report["Status"] === "OK" &&
        report["Result"] === expected &&
        AssociationQ[report["Summary"]] &&
        AssociationQ[report["OrderingSummary"]]
    ],
    True,
    TestID -> "dnfReduceDiagnosticReport: matches dnfReduceSystem on small case"
]

Test[
    Module[{s, sys, report},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        sys = makeSystem[s];
        report = solversTools`Private`dnfReduceDiagnosticReport[sys, "Timeout" -> 0];
        report["Status"] === "Timeout" &&
        report["Result"] === $TimedOut &&
        AssociationQ[report["TimeoutLocation"]]
    ],
    True,
    TestID -> "dnfReduceDiagnosticReport: timeout returns diagnostic report"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        solversTools`Private`solutionResultKind[result]
    ],
    "Rules",
    TestID -> "solutionResultKind: chain-2v dnf result is Rules"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        solversTools`Private`solutionResultKind[result]
    ],
    "Branched",
    TestID -> "solutionResultKind: chain-3v-2exit dnf result is Branched"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[7, {{1, 100.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        solversTools`Private`solutionResultKind[result]
    ],
    "Branched",
    TestID -> "solutionResultKind: example-7 dnf result is Branched"
]

Test[
    Module[{s, sys, result, kind},
        s = getExampleScenario[12, {{1, 100.0}}, {{4, 0.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        kind = solversTools`Private`solutionResultKind[result];
        MemberQ[{"Branched", "Parametric"}, kind]
    ],
    True,
    TestID -> "solutionResultKind: example-12 dnf result is specific residual kind"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "dnfReduceSystem: chain 2-exits solution is valid"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[8, {{1, 80.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "dnfReduceSystem: y-network solution is valid"
]

Test[
    Module[{s, sys, result, rules, residual},
        s = getExampleScenario[7, {{1, 100.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        rules = If[ListQ[result], result, Lookup[result, "Rules", {}]];
        residual = If[AssociationQ[result], Lookup[result, "Equations", True], True];
        And[
            AssociationQ[result],
            (j[2, 4] /. rules) === 100 - j[2, 3],
            (j[3, "auxExit3"] /. rules) === j[2, 3],
            (j[4, "auxExit4"] /. rules) === 100 - j[2, 3],
            !FreeQ[residual, j[2, 3] == 55],
            isValidSystemSolution[sys, result]
        ]
    ],
    True,
    TestID -> "dnfReduceSystem: example-7 preserves edge-flow branches"
]

Test[
    Module[{s, sys, result, rules, residual},
        s = getExampleScenario[7, {{1, 100.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        rules = If[ListQ[result], result, Lookup[result, "Rules", {}]];
        residual = If[AssociationQ[result], Lookup[result, "Equations", True], True];
        And[
            AssociationQ[result],
            (u[1, 2] /. rules) === -100 + u[2, 1],
            !FreeQ[residual, u[2, 3] == 0],
            !FreeQ[residual, u[2, 4] == 10],
            isValidSystemSolution[sys, result]
        ]
    ],
    True,
    TestID -> "dnfReduceSystem: example-7 preserves value branches"
]

Test[
    Module[{s, sys, result, rules, residual},
        s = gridScenario[{3}, {{1, 5}}, {{2, 0}, {3, 10}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        rules = If[ListQ[result], result, Lookup[result, "Rules", {}]];
        residual = If[AssociationQ[result], Lookup[result, "Equations", True], True];
        And[
            AssociationQ[result],
            (j[2, 3] /. rules) === 0,
            (j[3, 2] /. rules) === 0,
            (u[2, 3] /. rules) === u[2, 3],
            TrueQ[Simplify[Implies[residual, u[2, 3] <= 10]]],
            isValidSystemSolution[sys, result]
        ]
    ],
    True,
    TestID -> "dnfReduceSystem: zero-flow edge leaves unused value constrained"
]

Test[
    Module[{s, sys, result, residual},
        s = gridScenario[{3}, {{1, 5}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        residual = If[AssociationQ[result], Lookup[result, "Equations", True], True];
        And[
            AssociationQ[result],
            !FreeQ[residual, j[2, "auxExit2"] == 5],
            !FreeQ[residual, j[2, "auxExit2"] == 0],
            isValidSystemSolution[sys, result]
        ]
    ],
    True,
    TestID -> "dnfReduceSystem: relaxed edge equation preserves chain branches"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2];
        sys = makeSystem[s];
        result = Quiet[dnfReduceSystem[sys], dnfReduceSystem::noncritical];
        FailureQ[result] && result["Tag"] === "dnfReduceSystem"
    ],
    True,
    TestID -> "dnfReduceSystem: non-critical congestion systems fail"
]

(* --- optimizedDNFReduceSystem tests --- *)

Test[
    NameQ["solversTools`optimizedDNFReduceSystem"],
    True,
    TestID -> "optimizedDNFReduceSystem: public symbol exists"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = optimizedDNFReduceSystem[sys];
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "optimizedDNFReduceSystem: chain 2-exits solution is valid"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[8, {{1, 80.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = optimizedDNFReduceSystem[sys];
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "optimizedDNFReduceSystem: y-network solution is valid"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[7, {{1, 100.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = optimizedDNFReduceSystem[sys];
        AssociationQ[result] && isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "optimizedDNFReduceSystem: example-7 returns valid exact residual family"
]

Test[
    Module[{s, sys, result, rules, residual},
        s = getExampleScenario[12, {{1, 100.0}}, {{4, 0.0}}];
        sys = makeSystem[s];
        result = optimizedDNFReduceSystem[sys];
        rules = If[ListQ[result], result, Lookup[result, "Rules", {}]];
        residual = If[AssociationQ[result], Lookup[result, "Equations", True], True];
        AssociationQ[result] &&
        FreeQ[rules, _Real] &&
        FreeQ[residual, _Real] &&
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "optimizedDNFReduceSystem: example-12 remains valid and exact"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2];
        sys = makeSystem[s];
        result = Quiet[optimizedDNFReduceSystem[sys], optimizedDNFReduceSystem::noncritical];
        FailureQ[result] && result["Tag"] === "optimizedDNFReduceSystem"
    ],
    True,
    TestID -> "optimizedDNFReduceSystem: non-critical congestion systems fail"
]

(* --- activeSetReduceSystem tests --- *)

Test[
    NameQ["solversTools`activeSetReduceSystem"],
    True,
    TestID -> "activeSetReduceSystem: public symbol exists"
]

Test[
    StringContainsQ[
        ToString[DownValues[solversTools`activeSetReduceSystem], InputForm],
        "branchStateReduceFromBranches"
    ],
    True,
    TestID -> "activeSetReduceSystem: uses distinct active-set branch path"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = activeSetReduceSystem[sys];
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "activeSetReduceSystem: chain 2-exits solution is valid"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[8, {{1, 80.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = activeSetReduceSystem[sys];
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "activeSetReduceSystem: y-network solution is valid"
]

Test[
    Module[{s, sys, result, rules, residual},
        s = getExampleScenario[7, {{1, 100.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = activeSetReduceSystem[sys];
        rules = If[ListQ[result], result, Lookup[result, "Rules", {}]];
        residual = If[AssociationQ[result], Lookup[result, "Equations", True], True];
        AssociationQ[result] &&
        FreeQ[rules, _Real] &&
        FreeQ[residual, _Real] &&
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "activeSetReduceSystem: example-7 returns valid exact residual family"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[12, {{1, 100.0}}, {{4, 0.0}}];
        sys = makeSystem[s];
        result = activeSetReduceSystem[sys];
        AssociationQ[result] && isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "activeSetReduceSystem: example-12 remains valid"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2];
        sys = makeSystem[s];
        result = Quiet[activeSetReduceSystem[sys], activeSetReduceSystem::noncritical];
        FailureQ[result] && result["Tag"] === "activeSetReduceSystem"
    ],
    True,
    TestID -> "activeSetReduceSystem: non-critical congestion systems fail"
]

Test[
    Module[{s, sys},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        sys = makeSystem[s];
        !isValidSystemSolution[
            sys,
            <|"Rules" -> {}, "Equations" -> (u[1, 2] == 0 && u[1, 2] == 1)|>
        ]
    ],
    True,
    TestID -> "isValidSystemSolution: rejects inconsistent residual"
]

(* --- booleanReduceSystem tests --- *)

Test[
    NameQ["solversTools`booleanReduceSystem"],
    True,
    TestID -> "booleanReduceSystem: public symbol exists"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = booleanReduceSystem[sys];
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "booleanReduceSystem: chain 2-exits solution is valid"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[8, {{1, 80.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = Quiet[booleanReduceSystem[sys], booleanReduceSystem::multisol];
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "booleanReduceSystem: y-network solution is valid"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2];
        sys = makeSystem[s];
        result = Quiet[booleanReduceSystem[sys], booleanReduceSystem::noncritical];
        FailureQ[result] && result["Tag"] === "booleanReduceSystem"
    ],
    True,
    TestID -> "booleanReduceSystem: non-critical congestion systems fail"
]

(* --- findInstanceSystem tests --- *)

Test[
    NameQ["solversTools`findInstanceSystem"],
    True,
    TestID -> "findInstanceSystem: public symbol exists"
]

Test[
    Module[{data, s, sys, result},
        data = <|
            "Vertices" -> {1, 2, 3},
            "Adjacency" -> {
                {0, 1, 0},
                {1, 0, 1},
                {0, 1, 0}
            },
            "Entries" -> {{1, 10}},
            "Exits" -> {{3, 0}},
            "Switching" -> {}
        |>;
        s = makeScenario[<|"Model" -> data|>];
        sys = makeSystem[s];
        result = findInstanceSystem[sys];
        MatchQ[result, {__Rule}] && isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "findInstanceSystem: chain 1-exit returns valid rules"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = findInstanceSystem[sys];
        MatchQ[result, {__Rule}] && isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "findInstanceSystem: chain 2-exits solution is valid"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[8, {{1, 80.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = findInstanceSystem[sys];
        MatchQ[result, {__Rule}] && isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "findInstanceSystem: y-network solution is valid"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2];
        sys = makeSystem[s];
        result = Quiet[findInstanceSystem[sys], findInstanceSystem::noncritical];
        FailureQ[result] && result["Tag"] === "findInstanceSystem"
    ],
    True,
    TestID -> "findInstanceSystem: non-critical congestion systems fail"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[8, {{1, 80.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = findInstanceSystem[sys, "Timeout" -> 0];
        AssociationQ[result] &&
        KeyExistsQ[result, "Rules"] &&
        Lookup[result, "Equations", Missing["KeyAbsent", "Equations"]] === False
    ],
    True,
    TestID -> "findInstanceSystem: timeout returns rules plus false residual"
]
