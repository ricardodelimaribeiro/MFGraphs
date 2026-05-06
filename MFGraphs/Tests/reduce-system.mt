(* Tests for reduceSystem *)

Test[
    NameQ["solversTools`reduceSystem"],
    True,
    TestID -> "reduceSystem: public symbol exists"
]

Test[
    NameQ["systemTools`buildBoundaryData"] &&
    NameQ["systemTools`buildFlowData"] &&
    NameQ["systemTools`buildComplementarityData"] &&
    NameQ["systemTools`buildHamiltonianData"],
    True,
    TestID -> "systemTools: modular system builders are public symbols"
]

Test[
    Module[{s, sys, data},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        sys = makeSystem[s];
        data = systemDataFlatten[sys];
        AssociationQ[data] &&
        KeyExistsQ[data, "EqGeneral"] &&
        KeyExistsQ[data, "IneqExitValues"] &&
        KeyExistsQ[data, "AltExitCond"] &&
        KeyExistsQ[data, "AltOptCond"] &&
        !KeyExistsQ[data, "BoundaryData"] &&
        !KeyExistsQ[data, "FlowData"]
    ],
    True,
    TestID -> "systemDataFlatten: exposes nested solver keys without typed sub-record wrappers"
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
        exitVals = Cases[systemData[sys, "IneqExitValues"], HoldPattern[_ <= v_] :> v];
        FreeQ[Join[entryVals, exitVals], _Real]
    ],
    True,
    TestID -> "reduceSystem: boundary rules are exactified (no Real coefficients)"
]

Test[
    Module[{s, sys, exitIneqs, altExit},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        exitIneqs = systemData[sys, "IneqExitValues"];
        altExit   = systemData[sys, "AltExitCond"];
        MemberQ[exitIneqs, u["auxExit2", 2] <= 0] &&
        MemberQ[exitIneqs, u["auxExit3", 3] <= 10] &&
        !MemberQ[exitIneqs, u[2, "auxExit2"] <= _] &&
        MemberQ[altExit, (j[2, "auxExit2"] == 0) || (u["auxExit2", 2] == 0)] &&
        MemberQ[altExit, (j[3, "auxExit3"] == 0) || (u["auxExit3", 3] == 10)]
    ],
    True,
    TestID -> "reduceSystem: IneqExitValues and AltExitCond use u[auxExit, vertex] orientation"
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
        result = Quiet[reduceSystem[sys], MFGraphs::noncritical];
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
    <|"Rules" -> {}, "Residual" -> False|>,
    TestID -> "branchStateReduceResult: rejects direct conflicting equalities"
]

Test[
    solversTools`Private`branchStateReduceResult[(x == 1 || x == 2) && x == 3, {x}],
    <|"Rules" -> {}, "Residual" -> False|>,
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
        !FreeQ[Lookup[result, "Residual", True], x + y == 1]
    ],
    True,
    TestID -> "branchStateReduceResult: symbolic RHS remains residual"
]

Test[
    solversTools`Private`solutionResultKind[
        solversTools`Private`branchStateReduceResult[x + y == 1, {x, y}]
    ],
    "ResidualLogic",
    TestID -> "solutionResultKind: symbolic RHS branch-state result is residual"
]

Test[
    solversTools`Private`branchStateReduceResult[x == y && y == 3, {x, y}],
    {x -> 3, y -> 3},
    TestID -> "branchStateReduceResult: chained equalities become ground rules"
]

Test[
    utilities`mergeRules[{x -> 1}, {x -> 2}],
    {x -> 2},
    TestID -> "mergeRules: duplicate lhs keeps latest rule"
]

Test[
    solversTools`Private`solutionResultKind[{j[1, 2] -> 10}],
    "Rules",
    TestID -> "solutionResultKind: rule list is Rules"
]

Test[
    solversTools`Private`solutionResultKind[<|"Rules" -> {j[1, 2] -> 10}, "Residual" -> True|>],
    "Rules",
    TestID -> "solutionResultKind: true residual association is Rules"
]

Test[
    solversTools`Private`solutionResultKind[<|"Rules" -> {}, "Residual" -> False|>],
    "NoSolution",
    TestID -> "solutionResultKind: false residual association is NoSolution"
]

Test[
    solversTools`Private`solutionResultKind[<|"Rules" -> {}, "Residual" -> (j[1, 2] == 0 || j[1, 2] == 1)|>],
    "Branched",
    TestID -> "solutionResultKind: residual with Or is Branched"
]

Test[
    solversTools`Private`solutionResultKind[<|"Rules" -> {}, "Residual" -> j[1, 2, 3] >= 0|>],
    "Rules",
    TestID -> "solutionResultKind: transition-only residual is Rules"
]

Test[
    solversTools`Private`solutionResultKind[<|"Rules" -> {}, "Residual" -> (j[1, 2, 3] == 0 || j[1, 2, 3] == 1)|>],
    "Rules",
    TestID -> "solutionResultKind: transition-only Or is not Branched"
]

Test[
    solversTools`Private`solutionResultKind[<|"Rules" -> {}, "Residual" -> (j[1, 2] >= 0 && u[1, 2] <= 1)|>],
    "Underdetermined",
    TestID -> "solutionResultKind: tracked residual without Or is Underdetermined"
]

Test[
    solversTools`Private`solutionResultKind[
        <|"Rules" -> {}, "Residual" -> (j[1, 2] >= 0 && (j[1, 2, 3] == 0 || j[1, 2, 3] == 1))|>
    ],
    "Underdetermined",
    TestID -> "solutionResultKind: mixed residual classifies by primary variables"
]

Test[
    solversTools`Private`solutionResultKind[
        <|"Rules" -> {}, "Residual" -> ((j[1, 2] == 0 && j[1, 2, 3] >= 0) || j[1, 2] == 1)|>
    ],
    "Branched",
    TestID -> "solutionResultKind: primary Or remains Branched"
]

Test[
    Module[{diag},
        diag = solversTools`Private`dnfResidualDiagnostics[
            <|"Rules" -> {}, "Residual" -> ((j[1, 2] == 0 && (u[1, 2] == 0 || u[1, 2] == 1)) || j[1, 2] == 1)|>
        ];
        diag["TopLevelBranchCount"] === 2 &&
        diag["NestedOrQ"] === True &&
        diag["TrackedVariablesQ"] === True &&
        diag["PrimaryVariables"] === {j[1, 2], u[1, 2]} &&
        diag["TransitionFlowVariables"] === {}
    ],
    True,
    TestID -> "dnfResidualDiagnostics: reports branches nested Or and tracked variables"
]

Test[
    Module[{diag},
        diag = solversTools`Private`dnfResidualDiagnostics[
            <|"Rules" -> {}, "Residual" -> (j[1, 2] >= 0 && j[1, 2, 3] >= 0)|>
        ];
        diag["PrimaryVariables"] === {j[1, 2]} &&
        diag["TransitionFlowVariables"] === {j[1, 2, 3]}
    ],
    True,
    TestID -> "dnfResidualDiagnostics: separates primary and transition variables"
]

Test[
    Module[{s, sys, jts, diag},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        jts = systemData[sys, "Jts"];
        diag = solversTools`Private`solutionVariableDiagnostics[
            sys,
            Thread[jts -> Range[Length[jts]]]
        ];
        diag["PrimaryResultKind"] === "Rules" &&
        diag["TransitionFlowStatus"] === "Unique" &&
        diag["TransitionFlowCount"] === Length[jts] &&
        diag["TransitionFlowRuleCount"] === Length[jts] &&
        diag["TransitionFlowResidualCount"] === 0 &&
        diag["ResidualTransitionFlows"] === {}
    ],
    True,
    TestID -> "solutionVariableDiagnostics: transition flow rules are unique"
]

Test[
    Module[{s, sys, jt, diag},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}];
        sys = makeSystem[s];
        jt = First[systemData[sys, "Jts"]];
        diag = solversTools`Private`solutionVariableDiagnostics[
            sys,
            <|"Rules" -> {}, "Residual" -> jt >= 0|>
        ];
        diag["PrimaryResultKind"] === "Rules" &&
        diag["TransitionFlowStatus"] === "Underdetermined" &&
        diag["TransitionFlowResidualCount"] === 1 &&
        diag["ResidualTransitionFlows"] === {jt}
    ],
    True,
    TestID -> "solutionVariableDiagnostics: residual transition flow is underdetermined"
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
    (* u[auxExit3,3] is free in [0,10]: no flow reaches exit 3,
       so the exit value is unconstrained by complementarity. The rest of the
       solution is unaffected by its value — genuinely underdetermined. *)
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        solversTools`Private`solutionResultKind[result]
    ],
    "Underdetermined",
    TestID -> "solutionResultKind: chain-3v-2exit dnf result is Underdetermined"
]

Test[
    Module[{s, sys, result},
        s = getExampleScenario[7, {{1, 100.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        solversTools`Private`solutionResultKind[result]
    ],
    "Rules",
    TestID -> "solutionResultKind: example-7 dnf result is Rules"
]

Test[
    (* EqGeneral is now unconditional: zero-flow edges force u[a,b]=u[b,a], which
       propagates through switching inequalities to pin values at unused exits.
       Example-12 is now fully determined ("Rules") rather than "Branched"/"Parametric". *)
    Module[{s, sys, result, kind},
        s = getExampleScenario[12, {{1, 100.0}}, {{4, 0.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        kind = solversTools`Private`solutionResultKind[result];
        kind === "Rules" && isValidSystemSolution[sys, result]
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
    Module[{s, sys, result, rules},
        s = getExampleScenario[7, {{1, 100.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        rules = If[ListQ[result], result, Lookup[result, "Rules", {}]];
        And[
            ListQ[rules],
            (j[2, 3] /. rules) === 55,
            (j[2, 4] /. rules) === 45,
            (j[3, "auxExit3"] /. rules) === 55,
            (j[4, "auxExit4"] /. rules) === 45,
            isValidSystemSolution[sys, result]
        ]
    ],
    True,
    TestID -> "dnfReduceSystem: example-7 fully determines edge flows"
]

Test[
    Module[{s, sys, result, rules},
        s = getExampleScenario[7, {{1, 100.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        rules = If[ListQ[result], result, Lookup[result, "Rules", {}]];
        And[
            ListQ[rules],
            (u[1, 2] /. rules) === 55,
            (u[2, 3] /. rules) === 0,
            (u[2, 4] /. rules) === 10,
            isValidSystemSolution[sys, result]
        ]
    ],
    True,
    TestID -> "dnfReduceSystem: example-7 fully determines value functions"
]

Test[
    Module[{s, sys, result, rules},
        s = gridScenario[{3}, {{1, 5}}, {{2, 0}, {3, 10}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        rules = If[ListQ[result], result, Lookup[result, "Rules", {}]];
        And[
            ListQ[rules],
            (j[2, 3] /. rules) === 0,
            (j[3, 2] /. rules) === 0,
            (u[2, 3] /. rules) === 0,
            isValidSystemSolution[sys, result]
        ]
    ],
    True,
    TestID -> "dnfReduceSystem: cheaper-exit chain fully determines zero flow and value"
]

Test[
    Module[{s, sys, result, rules},
        s = gridScenario[{3}, {{1, 5}}, {{2, 10}, {3, 0}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        rules = If[ListQ[result], result, Lookup[result, "Rules", {}]];
        And[
            ListQ[rules],
            (j[2, 3] /. rules) === 5,
            (j[2, "auxExit2"] /. rules) === 0,
            isValidSystemSolution[sys, result]
        ]
    ],
    True,
    TestID -> "dnfReduceSystem: cheaper-remote-exit chain fully determines flow"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2];
        sys = makeSystem[s];
        result = Quiet[dnfReduceSystem[sys], MFGraphs::noncritical];
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
    (* With unconditional EqGeneral the solution is now a fully-determined List of
       rules rather than an Association with residual branches. *)
    Module[{s, sys, result},
        s = getExampleScenario[7, {{1, 100.0}}, {{3, 0.0}, {4, 10.0}}];
        sys = makeSystem[s];
        result = optimizedDNFReduceSystem[sys];
        !FailureQ[result] && isValidSystemSolution[sys, result]
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
        residual = If[AssociationQ[result], Lookup[result, "Residual", True], True];
        !FailureQ[result] &&
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
        result = Quiet[optimizedDNFReduceSystem[sys], MFGraphs::noncritical];
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
        residual = If[AssociationQ[result], Lookup[result, "Residual", True], True];
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
        !FailureQ[result] && isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "activeSetReduceSystem: example-12 remains valid"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2];
        sys = makeSystem[s];
        result = Quiet[activeSetReduceSystem[sys], MFGraphs::noncritical];
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
            <|"Rules" -> {}, "Residual" -> (u[1, 2] == 0 && u[1, 2] == 1)|>
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
        result = Quiet[booleanReduceSystem[sys], MFGraphs::noncritical];
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
        result = Quiet[findInstanceSystem[sys], MFGraphs::noncritical];
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
        Lookup[result, "Residual", Missing["KeyAbsent", "Residual"]] === False
    ],
    True,
    TestID -> "findInstanceSystem: timeout returns rules plus false residual"
]

(* --- ZeroSwitchUEqualities --- *)

Test[
    Module[{s, sys},
        s = gridScenario[{2, 3}, {{1, 100}}, {{6, 0}}];
        sys = makeSystem[s];
        Length[systemData[sys, "ZeroSwitchUEqualities"]] > 0
    ],
    True,
    TestID -> "ZeroSwitchUEqualities: grid-2x3 produces nontrivial equality rules"
]

Test[
    Module[{s, sys},
        s = gridScenario[{2, 3}, {{1, 100}}, {{6, 0}}];
        sys = makeSystem[s];
        Length[systemData[sys, "ZeroSwitchUEqualities"]] == 8
    ],
    True,
    TestID -> "ZeroSwitchUEqualities: grid-2x3 has 8 equality rules (one per non-canonical u at each vertex)"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120}}, {{2, 0}, {3, 10}}];
        sys = makeSystem[s];
        result = dnfReduceSystem[sys];
        isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "ZeroSwitchUEqualities: chain-2-exits solution remains valid after substitution"
]

Test[
    Module[{sSC, sNoSC, sysSC, sysNoSC},
        sSC   = gridScenario[{3}, {{1, 120}}, {{3, 0}}, {{1, 2, 3, 5}}];
        sNoSC = gridScenario[{3}, {{1, 120}}, {{3, 0}}];
        sysSC   = makeSystem[sSC];
        sysNoSC = makeSystem[sNoSC];
        Length[systemData[sysSC, "ZeroSwitchUEqualities"]] <
        Length[systemData[sysNoSC, "ZeroSwitchUEqualities"]]
    ],
    True,
    TestID -> "ZeroSwitchUEqualities: nonzero switching cost suppresses equality rule"
]

Test[
    Module[{s, sys},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        sys = makeSystem[s];
        ListQ[systemData[sys, "ZeroSwitchUEqualities"]]
    ],
    True,
    TestID -> "ZeroSwitchUEqualities: key exists and is a list"
]

Test[
    Module[{s, sys, b, km, vars},
        s = gridScenario[{3}, {{1, 10}}, {{2, 0}, {3, 5}}];
        sys = makeSystem[s];
        {b, km, vars} = getKirchhoffLinearSystem[sys];
        VectorQ[vars, MatchQ[#, j[__]] &] &&
        MatrixQ[km] &&
        Dimensions[km][[2]] === Length[vars] &&
        Length[b] === Dimensions[km][[1]]
    ],
    True,
    TestID -> "getKirchhoffLinearSystem: nontrivial chain returns aligned matrix dimensions"
]

(* --- Additive legacy-inspired diagnostics and opt-in solvers --- *)

Test[
    NameQ["solversTools`solutionReport"] &&
    NameQ["solversTools`solutionBranchCostReport"] &&
    NameQ["solversTools`directCriticalSystem"] &&
    NameQ["solversTools`flowFirstCriticalSystem"],
    True,
    TestID -> "legacy imports: public additive helper symbols exist"
]

Test[
    Module[{s, sys, sol, before, report},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        sys = makeSystem[s];
        sol = dnfReduceSystem[sys];
        before = sol;
        report = solutionReport[sys, sol];
        sol === before &&
        AssociationQ[report] &&
        report["ResultKind"] === "Rules" &&
        TrueQ[Lookup[report["ValidationReport"], "Valid", False]] &&
        NumericQ[report["KirchhoffResidual"]] &&
        AssociationQ[report["TransitionFlowDiagnostics"]] &&
        ListQ[report["ResidualTransitionFlows"]]
    ],
    True,
    TestID -> "solutionReport: reports diagnostics without rewriting solution"
]

Test[
    Module[{s, sys, baseRules, sol, report, totals},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        sys = makeSystem[s];
        baseRules = dnfReduceSystem[sys];
        baseRules = If[ListQ[baseRules], baseRules, Lookup[baseRules, "Rules", {}]];
        baseRules = DeleteCases[baseRules, HoldPattern[j[1, 2] -> _]];
        sol = <|"Rules" -> baseRules, "Residual" -> (j[1, 2] == 5 || j[1, 2] == 10)|>;
        report = solutionBranchCostReport[sys, sol];
        totals = Lookup[#, "TotalObjective"] & /@ report["Branches"];
        report["BranchCount"] === 2 &&
        Length[report["BestBranches"]] >= 1 &&
        TrueQ[First[totals] <= Last[totals]]
    ],
    True,
    TestID -> "solutionBranchCostReport: ranks residual branches by objective"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        sys = makeSystem[s];
        result = directCriticalSystem[sys];
        AssociationQ[result] && isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "directCriticalSystem: equal-exit zero-switching critical case succeeds"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 10}}, {{2, 0}, {3, 5}}];
        sys = makeSystem[s];
        result = directCriticalSystem[sys];
        FailureQ[result] && result["Tag"] === "directCriticalSystem" &&
        result["Reason"] === "ExitCostsNotEqualNumeric"
    ],
    True,
    TestID -> "directCriticalSystem: unequal exits fail explicitly"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}, {{1, 2, 3, 5}}];
        sys = makeSystem[s];
        result = directCriticalSystem[sys];
        FailureQ[result] && result["Tag"] === "directCriticalSystem" &&
        result["Reason"] === "NonZeroOrNonNumericSwitchingCosts"
    ],
    True,
    TestID -> "directCriticalSystem: nonzero switching costs fail explicitly"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 10}}, {{3, 0}}, {{1, 2, 3, Infinity}}];
        sys = makeSystem[s];
        result = directCriticalSystem[sys];
        FailureQ[result] && result["Tag"] === "directCriticalSystem" &&
        result["Reason"] === "NonZeroOrNonNumericSwitchingCosts"
    ],
    True,
    TestID -> "directCriticalSystem: infinite switching costs fail explicit zero-cost precondition"
]

Test[
    Module[{s, sys, result},
        s = gridScenario[{3}, {{1, 120}}, {{2, 0}, {3, 10}}];
        sys = makeSystem[s];
        result = flowFirstCriticalSystem[sys];
        AssociationQ[result] && isValidSystemSolution[sys, result]
    ],
    True,
    TestID -> "flowFirstCriticalSystem: small critical unequal-exit case validates"
]

Test[
    Module[{s, sys, defaultResult, directResult},
        s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
        sys = makeSystem[s];
        defaultResult = solveScenario[s];
        directResult = directCriticalSystem[sys];
        defaultResult === dnfReduceSystem[sys] &&
        directResult =!= defaultResult
    ],
    True,
    TestID -> "opt-in solvers: solveScenario default remains dnfReduceSystem"
]
