Needs["MFGraphs`"];
Needs["validationTools`"];

Test[
    Module[{sys, rules, report},
        sys = makeSystem[gridScenario[{2}, {}, {}]];
        rules = {j[1, 2] -> 0, j[2, 1] -> 0, u[1, 2] -> I, u[2, 1] -> I};
        report = exactSolutionReport[sys, rules];
        Lookup[report, {"Nonempty", "Soundness", "ResultType"}]],
    {False, "Vacuous", "EmptyReturnedSet"},
    TestID -> "exact validation: cancelling complex values do not make a real solution"
]
Test[
    Module[{sys, rules, report},
        sys = makeSystem[gridScenario[{2}, {}, {}]];
        rules = {j[1, 2] -> 0, j[2, 1] -> 0, u[1, 2] -> 1. + I, u[2, 1] -> 1. + I};
        report = exactSolutionReport[sys, rules];
        Lookup[report, {"ExactOutput", "ResultType"}]],
    {False, "NumericalApproximation"},
    TestID -> "exact validation: inexact complex atoms retain numerical provenance"
]
Test[
    Module[{sys, rules, report, parameter},
        sys = makeSystem[gridScenario[{2}, {}, {}]];
        rules = {j[1, 2] -> 0, j[2, 1] -> 0,
            u[1, 2] -> Sqrt[parameter], u[2, 1] -> Sqrt[parameter]};
        report = exactSolutionReport[sys, rules, "CheckCompleteness" -> False];
        {report["Soundness"], Resolve[ForAll[parameter,
            Equivalent[report["ParameterDomain"], parameter >= 0]], Reals]}],
    {"Proved", True},
    TestID -> "exact validation: real image restricts external family parameter domain"
]
Test[
    Module[{sys, rules},
        sys = makeSystem[gridScenario[{2}, {}, {}]];
        rules = {j[1, 2] -> 0, j[2, 1] -> 0, u[1, 2] -> $Failed, u[2, 1] -> $Failed};
        exactSolutionReport[sys, rules]["ResultType"]],
    "Unresolved", TestID -> "exact validation: cancelling backend failure cannot become a solution"
]
Test[
    Module[{sys, report},
        sys = makeSystem[gridScenario[{2}, {}, {}]];
        report = exactSolutionReport[sys, <|"Rules" -> {j[1, 2] -> 0, j[2, 1] -> 0},
            "Residual" -> (u[1, 2] == u[2, 1] || u[1, 2] == u[2, 1] + 1)|>];
        report["Soundness"]],
    "Disproved", TestID -> "exact validation: every returned branch must satisfy the original system"
]
Test[
    Module[{sys},
        sys = makeSystem[gridScenario[{2}, {}, {}]];
        exactSolutionReport[sys, {u[1, 2] -> u[2, 1], u[2, 1] -> u[1, 2]}]["ResultType"]],
    "Unresolved", TestID -> "exact validation: cyclic replacement rules remain unresolved"
]
Test[
    Module[{sys, data},
        sys = makeSystem[gridScenario[{2}, {}, {}]];
        data = systemData[sys];
        exactSolutionReport[mfgSystem[Append[data, "Us" -> {}]], {}]["ResultType"]],
    "InvalidSystem", TestID -> "exact validation: original unknowns cannot disappear from variable accounting"
]
Test[
    Module[{sys},
        sys = makeSystem[makeScenario[<|"Model" -> <|
            "Graph" -> Graph[{1 <-> 2, 2 <-> 3, 2 <-> 4}],
            "Entries" -> {{1, 1}}, "Exits" -> {{3, 0}},
            "Switching" -> {{1, 2, 3, 1.0}, {1, 2, 4, 0}, {4, 2, 3, 0}}|>|>]];
        exactSolutionReport[sys, {}]["ExactInput"]],
    False, TestID -> "exact validation: metric closure cannot erase inexact source provenance"
]

Test[
    Resolve[ForAll[{forward, backward, net}, Equivalent[
        forward >= 0 && backward >= 0 && (forward == 0 || backward == 0) && forward - backward == net,
        forward == Max[net, 0] && backward == Max[-net, 0]]], Reals],
    True, TestID -> "net-flow reduction: identity for every real net current"
]
Test[
    Module[{sys, data, rules}, sys = makeSystem[gridScenario[{2}, {{1, 10}}, {{2, 0}}]];
        data = systemDataFlatten[sys]; rules = {j[1, 2] -> 10 + j[2, 1]};
        solversTools`Private`fixedNetFlowEqualities[Append[data, "AltFlows" -> {}], rules]],
    {}, TestID -> "net-flow reduction: missing complementarity cannot justify pinning"
]
Test[
    Module[{sys, data, rules}, sys = makeSystem[gridScenario[{2}, {{1, 10}}, {{2, 0}}]];
        data = systemDataFlatten[sys]; rules = {j[1, 2] -> 10. + j[2, 1]};
        solversTools`Private`fixedNetFlowEqualities[data, rules]],
    {}, TestID -> "net-flow reduction: approximate current is not an exact pin"
]
Do[
    With[{scen = pair[[2]], label = pair[[1]]}, Test[
        Module[{sys, result, report}, sys = makeSystem[scen];
            result = linearNetReduceSystem[sys];
            report = exactSolutionReport[sys, result, "TimeLimit" -> 5];
            Lookup[report, {"Soundness", "Completeness"}]],
        {"Proved", "Proved"}, TestID -> "net-flow reduction: original-system equivalence " <> label
    ]],
    {pair, {{"reverse-current", gridScenario[{2}, {{2, 7/3}}, {{1, 0}}]},
            {"two-exits", gridScenario[{3}, {{2, 80}}, {{1, 0}, {3, 10}}]},
            {"cycle", cycleScenario[3, {{1, 10}}, {{3, 0}}]}}}
];

Test[
    FailureQ[makeSystem[gridScenario[{3}, {{1, 10}}, {{3, 0}},
        {{1, 2, 3, Function[flow, flow^2]}}]]],
    True, TestID -> "system: nonlinear switching cannot enter linear elimination"
]
Test[
    solversTools`Private`rulesFromEqualities[{x^2 + y == 1}, {x, y}, {}],
    {}, TestID -> "elimination: do not drop nonlinear coefficient arrays"
]

Test[
    isValidSystemSolution[makeSystem[gridScenario[{2}, {{1, 10}}, {{2, 0}}]], {}],
    False, TestID -> "validation: empty rules cannot satisfy prescribed inflow"
]
Test[
    isValidSystemSolution[makeSystem[gridScenario[{2}, {{1, 10}}, {{2, 0}}]],
        {j["auxEntry1", 1] -> 10}],
    False, TestID -> "validation: unresolved blocks are not validation success"
]
Test[
    Module[{sys, sol, r}, sys = makeSystem[gridScenario[{2}, {{1, 10}}, {{2, 0}}]];
        sol = dnfReduceSystem[sys]; r = exactSolutionReport[sys, sol];
        Lookup[r, {"ResultType", "Soundness", "Completeness"}]],
    {"ExactDetermined", "Proved", "Proved"},
    TestID -> "exact validation: original-system equivalence for determined solution"
]
Test[
    Module[{sys, sol, r}, sys = makeSystem[gridScenario[{3}, {{1, 10}}, {{2, 0}, {3, 10}}]];
        sol = dnfReduceSystem[sys]; r = exactSolutionReport[sys, sol];
        Lookup[r, {"ResultType", "Soundness", "Completeness"}]],
    {"ExactParametricFamily", "Proved", "Proved"},
    TestID -> "exact validation: parametric exit value has proved domain"
]
Test[
    Module[{sys}, sys = makeSystem[gridScenario[{2}, {{1, 10}}, {{2, 0}}]];
        exactSolutionReport[sys, {}]["Soundness"]],
    "Disproved", TestID -> "exact validation: unconstrained output is not a solution"
]
Test[
    Module[{sys, sol}, sys = makeSystem[gridScenario[{2}, {{1, 10.}}, {{2, 0}}]];
        sol = dnfReduceSystem[sys]; exactSolutionReport[sys, sol]["ResultType"]],
    "NumericalApproximation", TestID -> "exact validation: rationalized boundary retains inexact provenance"
]
Test[
    Module[{sys}, sys = makeSystem[gridScenario[{2}, {{1, -1}}, {{2, 0}}]];
        exactSolutionReport[sys, <|"Rules" -> {}, "Residual" -> False|>]["ResultType"]],
    "ExactInfeasible", TestID -> "exact validation: infeasibility independently proved"
]
Test[
    Module[{sys}, sys = makeSystem[gridScenario[{2}, {{1, 10}}, {{2, 0}}]];
        exactSolutionReport[sys, <|"Rules" -> {}, "Residual" -> False|>]["Completeness"]],
    "Disproved", TestID -> "exact validation: false infeasibility claim rejected"
]
Test[
    Module[{sys}, sys = makeSystem[gridScenario[{2}, {{1, 10}}, {{2, 0}}]];
        exactSolutionReport[sys, $TimedOut]["ResultType"]],
    "Timeout", TestID -> "exact validation: timeout stays distinct from infeasibility"
]
Test[
    Module[{sys, sol, report, free, parameter, renamed},
        sys = makeSystem[gridScenario[{3}, {{1, 10}}, {{2, 0}, {3, 10}}]];
        sol = dnfReduceSystem[sys];
        report = exactSolutionReport[sys, sol];
        free = First[report["FreeVariables"]];
        renamed = <|"Rules" -> Append[sol["Rules"] /. free -> parameter, free -> parameter],
            "Residual" -> (sol["Residual"] /. free -> parameter)|>;
        Lookup[exactSolutionReport[sys, renamed], {"Soundness", "Completeness"}]],
    {"Proved", "Proved"}, TestID -> "exact validation: external family parameter is existential in completeness"
]
Test[
    Module[{sys}, sys = makeSystem[gridScenario[{2}, {{1, 10}}, {{2, 0}}, {}, 2]];
        exactSolutionReport[sys, {}]["ResultType"]],
    "UnsupportedModel", TestID -> "exact validation: noncritical model cannot receive critical certification"
]
Test[
    Module[{sys}, sys = makeSystem[gridScenario[{3}, {{2, 80}}, {{1, 0}, {3, 10}},
        {{1, 2, 3, 2}, {3, 2, 1, 2}}]];
        findInstanceSystem[sys, "Timeout" -> 0]],
    $TimedOut, TestID -> "findInstance: time limit is not an infeasibility proof"
]
