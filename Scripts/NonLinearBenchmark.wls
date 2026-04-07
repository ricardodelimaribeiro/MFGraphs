#!/usr/bin/env wolframscript
(* NonLinearBenchmark.wls — Comprehensive non-linear case testing (alpha ≠ 1) *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

Print["MFGraphs Non-Linear Benchmark Suite"];
Print["===================================="];
Print["Date: ", DateString[]];
Print["Mathematica: ", $Version];
Print[];

(* --- Setup --- *)
PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
MFGraphs`$MFGraphsVerbose = False;
Get[FileNameJoin[{Directory[], "Scripts", "BenchmarkHelpers.wls"}]];

(* --- Non-Linear Test Cases --- *)

$NonLinearTestCases = <|
    (* Test Suite 1: Small uniform exponents *)
    "Y7_alpha2" -> <|
        "description" -> "Y-Network, quadratic congestion (α=2)",
        "data" -> GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0},
        "alpha" -> Function[edge, 2],
        "expectedBehavior" -> "Flows more evenly distributed due to quadratic penalty"
    |>,
    "Y7_alpha125" -> <|
        "description" -> "Y-Network, sub-linear congestion (α=1.25)",
        "data" -> GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0},
        "alpha" -> Function[edge, 1.25],
        "expectedBehavior" -> "Intermediate between linear and moderate non-linearity"
    |>,
    "Y7_alpha15" -> <|
        "description" -> "Y-Network, moderate non-linearity (α=1.5)",
        "data" -> GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0},
        "alpha" -> Function[edge, 1.5],
        "expectedBehavior" -> "Intermediate between linear and quadratic"
    |>,

    (* Test Suite 2: Larger networks *)
    "Case12_alpha25" -> <|
        "description" -> "Attraction Network, high non-linearity (α=2.5)",
        "data" -> GetExampleData[12] /. {I1 -> 100, U1 -> 0},
        "alpha" -> Function[edge, 2.5],
        "expectedBehavior" -> "Strong discouragement of high flows"
    |>,
    "Braess_alpha15" -> <|
        "description" -> "Braess Paradox, moderate non-linearity (α=1.5)",
        "data" -> GetExampleData["Braess split"] /. {I1 -> 100, U1 -> 0},
        "alpha" -> Function[edge, 1.5],
        "expectedBehavior" -> "Non-linear effects on paradox equilibrium"
    |>,

    (* Test Suite 3: Edge-dependent exponents *)
    "Y7_heterogeneous" -> <|
        "description" -> "Y-Network with heterogeneous α by edge",
        "data" -> GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0},
        "alpha" -> Function[edge,
            Switch[edge,
                {1, 2}, 0.8,   (* entrance: low congestion *)
                {2, 3}, 2.0,   (* middle: severe *)
                {2, 4}, 1.5,   (* branch: moderate *)
                _, 1.0
            ]
        ],
        "expectedBehavior" -> "Flows avoid high-congestion edge {2,3}"
    |>,

    (* Test Suite 4: Extreme non-linearity *)
    "Case12_alpha3" -> <|
        "description" -> "Attraction Network, cubic costs (α=3)",
        "data" -> GetExampleData[12] /. {I1 -> 50, U1 -> 0},
        "alpha" -> Function[edge, 3],
        "expectedBehavior" -> "Very steep cost penalty, extreme flow avoidance"
    |>
|>;

Print["Test cases defined: ", Length[$NonLinearTestCases]];
Print[];

(* --- Solver Execution --- *)

SolveNonLinearCase[caseKey_String, spec_Association] :=
    Module[{data, d2e, d2eAlpha1, baselineAlpha1, nonlinear, monotone, caseResults},
        Print["  Case: ", spec["description"]];
        data = spec["data"];
        Print["    Converting to equations (for α=", ToString[spec["alpha"][{1,2}]], ")..."];
        d2e = DataToEquations[data];
        Print["    Equations ready."];

        (* Baseline: CriticalCongestionSolver with α=1 *)
        Print["    Running baseline: CriticalCongestionSolver (α=1)..."];
        baselineAlpha1 = Quiet @ Check[
            TimeConstrained[
                CriticalCongestionSolver[d2e],
                60,
                <|"Feasibility" -> "Timeout", "Message" -> "CriticalSolver timeout"|>
            ],
            <|"Feasibility" -> "Failed", "Message" -> "CriticalSolver failed"|>
        ];
        Print["    Baseline done."];

        (* NonLinearSolver with custom alpha *)
        Print["    Running NonLinearSolver (α≠1)..."];
        nonlinear = Quiet @ Check[
            TimeConstrained[
                NonLinearSolver[
                    d2e,
                    "CongestionExponentFunction" -> spec["alpha"],
                    "MaxIterations" -> 15,
                    "Tolerance" -> 10^-6
                ],
                60,
                <|"Feasibility" -> "Timeout", "Message" -> "NonLinearSolver timeout"|>
            ],
            <|"Feasibility" -> "Failed", "Message" -> "NonLinearSolver failed"|>
        ];
        Print["    NonLinear done."];

        (* MonotoneSolver with custom alpha *)
        Print["    Running MonotoneSolver (α≠1)..."];
        monotone = Quiet @ Check[
            TimeConstrained[
                MonotoneSolver[
                    d2e,
                    "CongestionExponentFunction" -> spec["alpha"],
                    "ResidualTolerance" -> 10^-5,
                    "MaxTime" -> 20,
                    "MaxSteps" -> 1500
                ],
                60,
                <|"Feasibility" -> "Timeout", "Message" -> "MonotoneSolver timeout"|>
            ],
            <|"Feasibility" -> "Failed", "Message" -> "MonotoneSolver failed"|>
        ];
        Print["    Monotone done."];

        (* Helper: Compare flow solutions from ComparableFlowVector *)
        FlowDifference[sol1_, sol2_] := Module[{vec1, vec2, diffs},
            vec1 = Lookup[sol1, "ComparableFlowVector", Missing["N/A"]];
            vec2 = Lookup[sol2, "ComparableFlowVector", Missing["N/A"]];
            If[And[ListQ[vec1], ListQ[vec2], Length[vec1] === Length[vec2]],
                diffs = Abs[vec1 - vec2];
                Max[diffs],
                Missing["N/A"]
            ]
        ];

        (* Extract comparison data: compare α=1 baseline with α≠1 solvers *)
        caseResults = <|
            "Case" -> caseKey,
            "Description" -> spec["description"],
            "Alpha1Residual" -> Lookup[baselineAlpha1, "NonLinearResidual",
                Lookup[baselineAlpha1, "KirchhoffResidual", Missing["N/A"]]],
            "Alpha1Feasible" -> ToString[IsFeasible[baselineAlpha1]],
            "NonLinearResidual" -> Lookup[nonlinear, "NonLinearResidual",
                Lookup[nonlinear, "KirchhoffResidual", Missing["N/A"]]],
            "MonotoneResidual" -> Lookup[monotone, "NonLinearResidual",
                Lookup[monotone, "KirchhoffResidual", Missing["N/A"]]],
            "NonLinearFeasible" -> ToString[IsFeasible[nonlinear]],
            "MonotoneFeasible" -> ToString[IsFeasible[monotone]],
            "FlowDiff_Alpha1_NL" -> FlowDifference[baselineAlpha1, nonlinear],
            "FlowDiff_NL_Monotone" -> FlowDifference[nonlinear, monotone],
            "NonLinearIterations" -> Lookup[Lookup[nonlinear, "Convergence", <||>], "Iterations", Missing["N/A"]],
            "MonotoneIterations" -> Lookup[Lookup[monotone, "Convergence", <||>], "Iterations", Missing["N/A"]]
        |>;

        caseResults
    ];

(* --- Run All Cases --- *)

Print["Running ", Length[$NonLinearTestCases], " non-linear benchmark cases..."];
Print[];

allResultsList = Normal[$NonLinearTestCases];
allResults = Table[
    caseKey = allResultsList[[i, 1]];
    Print["[", i, "/", Length[allResultsList], "] ", caseKey];
    SolveNonLinearCase[caseKey, $NonLinearTestCases[caseKey]],
    {i, 1, Length[allResultsList]}
];

Print[];
Print["All cases completed."];
Print[];

(* --- Summary Table --- *)

Print["=== RESULTS SUMMARY ==="];
Print[];

summaryTable = Grid[
    Prepend[
        allResults,
        {"Case", "α=1\nResidual", "NL\nResidual", "Monotone\nResidual",
         "Flow Diff\n(α=1 vs NL)", "Flow Diff\n(NL vs Monotone)", "NL Feasible"}
    ],
    Frame -> All,
    ItemSize -> {Automatic, 1.5}
];
Print[summaryTable];

Print[];
Print["=== SOLVER EQUIVALENCE ==="];
Print[];

equivalenceCheck = Table[
    Module[{a1res, nlres, mres, dAlpha1NL, dNLMono},
        a1res = result["Alpha1Residual"];
        nlres = result["NonLinearResidual"];
        mres = result["MonotoneResidual"];
        dAlpha1NL = result["FlowDiff_Alpha1_NL"];
        dNLMono = result["FlowDiff_NL_Monotone"];
        <|
            "Case" -> result["Case"],
            "BaselineResidual" -> a1res,
            "FlowChangeAlpha1toNL" -> dAlpha1NL,
            "NLMonotoneAgreement" -> dNLMono,
            "AlphaSensitivity" -> If[
                NumericQ[dAlpha1NL] && dAlpha1NL > 10^-4,
                "Significant",
                "Minimal"
            ]
        |>
    ],
    {result, allResults}
];

Print[TableForm[equivalenceCheck]];

Print[];
Print["=== CONVERGENCE ANALYSIS ==="];
Print[];

convergenceData = Table[
    <|
        "Case" -> result["Case"],
        "Baseline(α=1)Feasible" -> result["Alpha1Feasible"],
        "NL(α≠1)Converged" -> result["NonLinearFeasible"],
        "Monotone(α≠1)Converged" -> result["MonotoneFeasible"],
        "NLIterations" -> result["NonLinearIterations"],
        "MonotoneIterations" -> result["MonotoneIterations"],
        "FlowSensitivity" -> If[
            NumericQ[result["FlowDiff_Alpha1_NL"]],
            ToString[NumberForm[result["FlowDiff_Alpha1_NL"], 3]],
            "N/A"
        ]
    |>,
    {result, allResults}
];

Print[TableForm[convergenceData]];

(* --- Export Results --- *)

Print[];
Print["Exporting results..."];

timestamp = StringReplace[DateString["ISODateTime"], ":" -> "-"];
resultsDir = FileNameJoin[{Directory[], "Results"}];

(* CSV export *)
csvPath = FileNameJoin[{resultsDir, "nonlinear_benchmark_" <> timestamp <> ".csv"}];
Export[csvPath, allResults, "CSV"];
Print["Results exported to: ", csvPath];

(* JSON export *)
jsonPath = FileNameJoin[{resultsDir, "nonlinear_benchmark_" <> timestamp <> ".json"}];
Export[jsonPath, allResults, "JSON"];
Print["JSON exported to: ", jsonPath];

(* Markdown report *)
mdPath = FileNameJoin[{resultsDir, "nonlinear_benchmark_" <> timestamp <> ".md"}];
mdCaseDetails = StringJoin[Table[
    "### " <> result["Case"] <> "\n" <>
    "Description: " <> result["Description"] <> "\n" <>
    "- Baseline (α=1) Residual: " <> ToString[result["Alpha1Residual"]] <> "\n" <>
    "- NonLinear (α≠1) Residual: " <> ToString[result["NonLinearResidual"]] <> "\n" <>
    "- Monotone (α≠1) Residual: " <> ToString[result["MonotoneResidual"]] <> "\n" <>
    "- Baseline Feasible: " <> result["Alpha1Feasible"] <> "\n" <>
    "- NonLinear Feasible: " <> result["NonLinearFeasible"] <> "\n" <>
    "- Monotone Feasible: " <> result["MonotoneFeasible"] <> "\n" <>
    "- Flow Change (α=1 → α≠1): " <> ToString[result["FlowDiff_Alpha1_NL"]] <> "\n" <>
    "- Flow Difference (NL vs Monotone): " <> ToString[result["FlowDiff_NL_Monotone"]] <> "\n\n",
    {result, allResults}
]];

mdContent = StringJoin[
    "# Non-Linear Benchmark Report\n\n",
    "Date: ", DateString[], "\n",
    "Mathematica: ", $Version, "\n\n",
    "## Test Cases\n\n",
    mdCaseDetails,
    "## Summary\n\n",
    "Total cases: ", ToString[Length[allResults]], "\n\n"
];

Export[mdPath, mdContent, "Text"];
Print["Markdown report: ", mdPath];

Print[];
Print["Benchmark complete."];

(* :!CodeAnalysis::EndBlock:: *)
