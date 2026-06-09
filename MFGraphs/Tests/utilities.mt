(* Unit tests for utilities.wl — typed-object infrastructure, rule management,
   numeric exactification, and the critical-congestion solver guard.

   These are pure functions feeding the whole solver, so they are tested
   directly here rather than only indirectly through the pipeline. *)

(* ------------------------------------------------------------------ *)
(* Typed-object infrastructure: mfgTypedQ / mfgData                    *)
(* ------------------------------------------------------------------ *)

Test[
    NameQ["utilities`mfgTypedQ"] && NameQ["utilities`mfgData"],
    True,
    TestID -> "utilities: mfgTypedQ / mfgData public symbols exist"
]

Test[
    Module[{head},
        mfgTypedQ[head[<|"a" -> 1|>], head] &&
        !mfgTypedQ[head[1], head] &&
        !mfgTypedQ[5, head]
    ],
    True,
    TestID -> "mfgTypedQ: true only for head[Association]"
]

Test[
    Module[{head},
        mfgData[head[<|"a" -> 1, "b" -> 2|>]] === <|"a" -> 1, "b" -> 2|> &&
        mfgData[head[<|"a" -> 1|>], "a"] === 1
    ],
    True,
    TestID -> "mfgData: returns full association and per-key value"
]

Test[
    Module[{head},
        mfgData[head[<|"a" -> 1|>], "missing"] === Missing["KeyAbsent", "missing"]
    ],
    True,
    TestID -> "mfgData: absent key returns Missing[KeyAbsent, key]"
]

(* ------------------------------------------------------------------ *)
(* Rule management: mergeRules / normalizeRules / extractRules         *)
(* ------------------------------------------------------------------ *)

Test[
    Module[{p, q, r},
        mergeRules[{p -> 1, q -> 2}, {q -> 9, r -> 3}] === {p -> 1, q -> 9, r -> 3}
    ],
    True,
    TestID -> "mergeRules: newRules wins per LHS, order preserved"
]

Test[
    Module[{p, q},
        (* a -> b chases b -> 3 to ground via ReplaceRepeated *)
        normalizeRules[{p -> q, q -> 3}] === {p -> 3, q -> 3}
    ],
    True,
    TestID -> "normalizeRules: rewrites RHS through the full rule set"
]

Test[
    Module[{p, q},
        extractRules[{p -> 1, 5, q -> 2, "noise"}] === {p -> 1, q -> 2}
    ],
    True,
    TestID -> "extractRules: selects only rules from a mixed list"
]

Test[
    Module[{p},
        extractRules[<|"Rules" -> {p -> 1}|>] === {p -> 1} &&
        extractRules[<|"NoRulesKey" -> 1|>] === {} &&
        extractRules[42] === {}
    ],
    True,
    TestID -> "extractRules: association and non-list inputs"
]

(* ------------------------------------------------------------------ *)
(* Numerics: roundValues                                               *)
(* ------------------------------------------------------------------ *)

Test[
    roundValues[2.0] == 2 &&
    roundValues["passthrough"] === "passthrough",
    True,
    TestID -> "roundValues: rounds numbers, leaves non-numerics untouched"
]

Test[
    Module[{x},
        Head[roundValues[x -> 2.0]] === Rule &&
        Head[roundValues[{1.0, 2.0}]] === List &&
        Length[roundValues[{1.0, 2.0}]] === 2 &&
        AssociationQ[roundValues[<|"a" -> 2.0|>]]
    ],
    True,
    TestID -> "roundValues: recurses through rules, lists, and associations"
]

(* ------------------------------------------------------------------ *)
(* Boundary exactification: exactBoundaryValue / exactBoundaryValues   *)
(* ------------------------------------------------------------------ *)

Test[
    exactBoundaryValue[3] === 3 &&
    exactBoundaryValue[1/2] === 1/2,
    True,
    TestID -> "exactBoundaryValue: exact numbers pass through unchanged"
]

Test[
    exactBoundaryValue[0.5] === 1/2 &&
    exactBoundaryValue[2.0] === 2,
    True,
    TestID -> "exactBoundaryValue: real machine numbers are rationalized exactly"
]

Test[
    FailureQ[exactBoundaryValue[1.0 + 2.0 I]],
    True,
    TestID -> "exactBoundaryValue: inexact complex value is a Failure"
]

Test[
    Module[{sym},
        FailureQ[exactBoundaryValue[sym]]
    ],
    True,
    TestID -> "exactBoundaryValue: non-numeric value is a Failure"
]

Test[
    exactBoundaryValues[{1, 0.5, 2.0}] === {1, 1/2, 2},
    True,
    TestID -> "exactBoundaryValues: maps over a list of valid values"
]

Test[
    Module[{sym},
        FailureQ[exactBoundaryValues[{1, sym, 0.5}]]
    ],
    True,
    TestID -> "exactBoundaryValues: returns the first Failure in the list"
]

(* ------------------------------------------------------------------ *)
(* Critical-congestion guard: criticalCongestionSystemQ /              *)
(* withCriticalCongestionGuard                                         *)
(* ------------------------------------------------------------------ *)

Test[
    Module[{s, sys},
        s   = gridScenario[{2}, {{1, 100}}, {{2, 0}}];
        sys = makeSystem[s];
        criticalCongestionSystemQ[sys]
    ],
    True,
    TestID -> "criticalCongestionSystemQ: default (Alpha == 1) system is critical"
]

Test[
    criticalCongestionSystemQ[mfgSystem[<||>]] === False,
    True,
    TestID -> "criticalCongestionSystemQ: system without HalfPairs is not critical"
]

Test[
    Module[{s, sys},
        s   = gridScenario[{2}, {{1, 100}}, {{2, 0}}];
        sys = makeSystem[s];
        (* HoldRest: the body is only evaluated on a critical system. *)
        withCriticalCongestionGuard[sys, "probe", 42] === 42
    ],
    True,
    TestID -> "withCriticalCongestionGuard: evaluates body on a critical system"
]

Test[
    FailureQ[Quiet[withCriticalCongestionGuard[mfgSystem[<||>], "probe", 42]]],
    True,
    TestID -> "withCriticalCongestionGuard: returns a Failure on a non-critical system"
]
