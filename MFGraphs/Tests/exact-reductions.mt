(* Logical equivalence and premise protection for exact flow propagation. *)
Needs["MFGraphs`"];

(* Raw balance metadata is deliberately absent: only enforced equations may
   justify additional constraints. Each test supplies its own real unknowns. *)
exactReductionBalanceData[equation_, nonnegative_List] := <|
    "EqBalanceSplittingFlows" -> equation,
    "EqBalanceGatheringFlows" -> True,
    "BalanceSplittingFlows" -> {},
    "BalanceGatheringFlows" -> {},
    "IneqJs" -> {},
    "IneqJts" -> nonnegative
|>;

Test[
    Module[{x, y, equation, signs, forced},
        equation = x/2 + 3 y/5 == 0;
        signs = {x >= 0, y >= 0};
        forced = solversTools`Private`zeroSumFlowEqualities[
            exactReductionBalanceData[equation, signs], {}, {}];
        Resolve[ForAll[{x, y}, Equivalent[
            And @@ signs && equation,
            And @@ signs && And @@ forced]], Reals]],
    True, TestID -> "zero-sum reduction: positive rational coefficients give an equivalent zero family"
]

Test[
    Module[{x, y, equation, signs, forced, data},
        equation = -2 x - 7 y == 0;
        signs = {x >= 0, y >= 0};
        data = Join[exactReductionBalanceData[True, signs],
            <|"EqBalanceGatheringFlows" -> equation|>];
        forced = solversTools`Private`zeroSumFlowEqualities[data, {}, {}];
        Resolve[ForAll[{x, y}, Equivalent[
            And @@ signs && equation,
            And @@ signs && And @@ forced]], Reals]],
    True, TestID -> "zero-sum reduction: negative coefficients work in gathering equations"
]

Test[
    Module[{x, y, forced},
        forced = solversTools`Private`zeroSumFlowEqualities[
            exactReductionBalanceData[2 x + 0 y == 0, {x >= 0, y >= 0}], {}, {}];
        (* y has zero coefficient and remains free, even though x is forced. *)
        Resolve[ForAll[{x, y}, Equivalent[And @@ forced, x == 0]], Reals]],
    True, TestID -> "zero-sum reduction: a zero coefficient does not force its variable"
]

Test[
    Module[{x, y},
        solversTools`Private`zeroSumFlowEqualities[
            exactReductionBalanceData[0 x + 0 y == 0, {x >= 0, y >= 0}], {}, {}]],
    {}, TestID -> "zero-sum reduction: an identically zero balance forces nothing"
]

Test[
    Module[{x, y},
        solversTools`Private`zeroSumFlowEqualities[
            exactReductionBalanceData[x - y == 0, {x >= 0, y >= 0}], {}, {}]],
    {}, TestID -> "zero-sum reduction: mixed signs retain nonzero equal flows"
]

Test[
    Module[{x, y},
        solversTools`Private`zeroSumFlowEqualities[
            exactReductionBalanceData[x + y == 1, {x >= 0, y >= 0}], {}, {}]],
    {}, TestID -> "zero-sum reduction: nonzero constant does not imply zero flows"
]

Test[
    Module[{x, y},
        solversTools`Private`zeroSumFlowEqualities[
            exactReductionBalanceData[x + y == 0, {x >= 0}], {}, {}]],
    {}, TestID -> "zero-sum reduction: every variable needs its nonnegativity premise"
]

Test[
    Module[{x, y, data},
        data = Join[exactReductionBalanceData[True, {x >= 0, y >= 0}], <|
            "BalanceSplittingFlows" -> {x},
            "BalanceGatheringFlows" -> {y}|>];
        solversTools`Private`zeroSumFlowEqualities[data, {}, {}]],
    {}, TestID -> "zero-sum reduction: metadata-only balances cannot create equations"
]

Do[
    With[{useRules = mode, label = If[mode, "accumulated rule", "proved equation"]}, Test[
        Module[{flow, x, y, equation, signs, forced},
            equation = flow == x + y;
            signs = {flow >= 0, x >= 0, y >= 0};
            forced = solversTools`Private`zeroSumFlowEqualities[
                exactReductionBalanceData[equation, signs],
                If[useRules, {}, {flow == 0}],
                If[useRules, {flow -> 0}, {}]];
            Resolve[ForAll[{flow, x, y}, Equivalent[
                And @@ signs && equation && flow == 0,
                And @@ signs && And @@ forced && flow == 0]], Reals]],
        True, TestID -> "zero-sum reduction: zero " <> label <> " preserves incident conservation"
    ]],
    {mode, {False, True}}
];

(* Compare the whole original formula to the reconstructed preprocessing output,
   without using a package solver or the exactSolutionReport reference path.
   Scalar variables keep quantification independent of compound j/u syntax. *)
Do[
    With[{scenarioSpec = pair[[2]], label = pair[[1]]}, Test[
        Module[{sys, names, blocks, original, inputs, reconstructed,
                vars, scalar, substitutions, proposition},
            sys = makeSystem[scenarioSpec];
            names = {"EqEntryIn", "EqBalanceSplittingFlows", "EqBalanceGatheringFlows",
                "EqGeneral", "IneqJs", "IneqJts", "IneqSwitchingByVertex",
                "IneqExitValues", "AltFlows", "AltTransitionFlows", "AltOptCond", "AltExitCond"};
            blocks = systemData[sys, #] & /@ names;
            original = And @@ (If[ListQ[#], And @@ #, #] & /@ blocks);
            inputs = solversTools`Private`buildFixedNetSolverInputs[sys];
            reconstructed = inputs[[1]] && And @@ (Equal @@@ inputs[[3]]);
            vars = DeleteDuplicates @ Join[systemData[sys, "Us"],
                systemData[sys, "Js"], systemData[sys, "Jts"]];
            scalar = Table[Unique["exactReductionVar$"], {Length[vars]}];
            substitutions = Thread[vars -> scalar];
            proposition = Equivalent[original, reconstructed] /. substitutions;
            With[{quantified = scalar, expression = proposition},
                TimeConstrained[Resolve[ForAll[quantified, expression], Reals], 5, $TimedOut]]],
        True, TestID -> "net preprocessing: original formula and reconstructed family are equivalent " <> label
    ]],
    {pair, {
        {"positive current", gridScenario[{2}, {{1, 7/3}}, {{2, 0}}]},
        {"negative current", gridScenario[{2}, {{2, 7/3}}, {{1, 0}}]},
        {"zero current and free values", gridScenario[{2}, {{1, 0}}, {{2, 0}}]}
    }}
];
