#!/usr/bin/env wolframscript
(* ::Package:: *)

(* VerifyNonLinearSolutions.wls

   Verifies whether solutions returned by NonLinearSolver and MonotoneSolver
   actually satisfy the non-linear MFG system equations.

   The KirchhoffResidual reported by both solvers only checks flow balance
   (a linear constraint). This script computes the actual non-linear residual
   Max[Abs[Nlhs - Nrhs]] against the solution to verify true correctness.
*)

SetDirectory[DirectoryName[$ScriptCommandLine[[1]]] <> "/.."];
Get["MFGraphs/MFGraphs.wl"];

Print["Non-Linear Solution Verification"];
Print["================================"];
Print[];

(* Define test cases - same as NonLinearBenchmark.wls *)
$TestCases = <|
    "Y7_alpha2" -> <|
        "data" -> GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0},
        "alpha" -> Function[edge, 2]
    |>,
    "Y7_alpha125" -> <|
        "data" -> GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0},
        "alpha" -> Function[edge, 1.25]
    |>,
    "Y7_alpha15" -> <|
        "data" -> GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0},
        "alpha" -> Function[edge, 1.5]
    |>,
    "Case12_alpha25" -> <|
        "data" -> GetExampleData[12] /. {I1 -> 50, U1 -> 0},
        "alpha" -> Function[edge, 2.5]
    |>,
    "Braess_alpha15" -> <|
        "data" -> GetExampleData["Braess Paradox"],
        "alpha" -> Function[edge, 1.5]
    |>,
    "Case12_alpha3" -> <|
        "data" -> GetExampleData[12] /. {I1 -> 50, U1 -> 0},
        "alpha" -> Function[edge, 3]
    |>
|>;

(* Helper: Extract solver solution association from result, trying multiple keys *)
ExtractSolution[result_] := Module[{sol},
    If[!AssociationQ[result], Return[Missing["NoResult"]]];
    sol = Lookup[result, "AssoNonCritical", Missing["NotFound"]];
    If[AssociationQ[sol] && Length[sol] > 0, Return[sol]];
    sol = Lookup[result, "AssoMonotone", Missing["NotFound"]];
    If[AssociationQ[sol] && Length[sol] > 0, Return[sol]];
    sol = Lookup[result, "Solution", Missing["NotFound"]];
    If[MatchQ[sol, {_Rule ..}], Return[Association @ sol]];
    If[AssociationQ[sol] && Length[sol] > 0, Return[sol]];
    Missing["NoSolution"]
];

(* Helper: Compute true non-linear residual against a solution *)
ComputeNonLinearResidual[d2e_, solution_] := Module[{Nlhs, Nrhs, lhsVals, rhsVals, diffs},
    If[!AssociationQ[solution], Return[Missing["NoSolution"]]];
    Nlhs = Lookup[d2e, "Nlhs", $Failed];
    Nrhs = Lookup[d2e, "Nrhs", $Failed];
    If[Nlhs === $Failed || Nrhs === $Failed, Return[Missing["NoNlhsRhs"]]];
    lhsVals = Quiet @ Check[N[Nlhs /. Normal[solution]], $Failed];
    rhsVals = Quiet @ Check[N[Nrhs /. Normal[solution]], $Failed];
    If[lhsVals === $Failed || rhsVals === $Failed, Return[Missing["EvalFailed"]]];
    diffs = Flatten[{lhsVals - rhsVals}];
    diffs = Select[diffs, NumericQ];
    If[diffs === {}, Return[Missing["NoNumeric"]]];
    Max[Abs[diffs]]
];

(* Helper: Compute Kirchhoff residual against a solution *)
ComputeKirchhoffResidual[d2e_, solution_] := Module[{B, KM, jj, jjVals, flowAssoc},
    If[!AssociationQ[solution], Return[Missing["NoSolution"]]];
    flowAssoc = Association @ KeySelect[solution, MatchQ[#, _MFGraphs`Private`j] &];
    {B, KM, jj} = MFGraphs`GetKirchhoffLinearSystem[d2e];
    jjVals = Lookup[flowAssoc, jj, Missing["NotInSolution"]];
    If[!VectorQ[jjVals, NumericQ], Return[Missing["NotNumeric"]]];
    N @ Norm[KM . jjVals - B, Infinity]
];

(* Run verification for each test case *)
allVerifications = {};

Do[
    Module[{spec, data, d2e, nlResult, monoResult, nlSol, monoSol,
            nlKirchhoff, nlNonLinear, monoKirchhoff, monoNonLinear,
            verification},

        Print["=== ", caseKey, " ==="];
        spec = $TestCases[caseKey];
        data = spec["data"];
        d2e = DataToEquations[data];

        (* Run NonLinearSolver with custom alpha *)
        Print["  Running NonLinearSolver..."];
        nlResult = Quiet @ TimeConstrained[
            NonLinearSolver[d2e, "CongestionExponentFunction" -> spec["alpha"],
                "MaxIterations" -> 15, "Tolerance" -> 10^-6],
            60,
            <|"Status" -> "Timeout"|>
        ];

        (* Run MonotoneSolver with custom alpha *)
        Print["  Running MonotoneSolver..."];
        monoResult = Quiet @ TimeConstrained[
            MonotoneSolver[d2e, "CongestionExponentFunction" -> spec["alpha"],
                "ResidualTolerance" -> 10^-5, "MaxTime" -> 20, "MaxSteps" -> 1500],
            60,
            <|"Status" -> "Timeout"|>
        ];

        (* Extract solutions (try AssoNonCritical, AssoMonotone, then Solution) *)
        nlSol = ExtractSolution[nlResult];
        monoSol = ExtractSolution[monoResult];

        Print["  NL solution keys: ", If[AssociationQ[nlSol], Length[nlSol], nlSol]];
        Print["  Mono solution keys: ", If[AssociationQ[monoSol], Length[monoSol], monoSol]];

        (* Compute residuals against the actual non-linear system *)
        Print["  Computing residuals..."];
        nlKirchhoff = ComputeKirchhoffResidual[d2e, nlSol];
        nlNonLinear = ComputeNonLinearResidual[d2e, nlSol];
        monoKirchhoff = ComputeKirchhoffResidual[d2e, monoSol];
        monoNonLinear = ComputeNonLinearResidual[d2e, monoSol];

        verification = <|
            "Case" -> caseKey,
            "Alpha" -> spec["alpha"][{1,2}],
            "NL_Kirchhoff" -> nlKirchhoff,
            "NL_NonLinear" -> nlNonLinear,
            "Mono_Kirchhoff" -> monoKirchhoff,
            "Mono_NonLinear" -> monoNonLinear,
            "NL_Solves_NonLinear" -> If[NumericQ[nlNonLinear], nlNonLinear < 10^-4, False],
            "Mono_Solves_NonLinear" -> If[NumericQ[monoNonLinear], monoNonLinear < 10^-4, False]
        |>;

        AppendTo[allVerifications, verification];

        Print["  NonLinear: Kirchhoff=", nlKirchhoff, " | NonLinearResidual=", nlNonLinear];
        Print["  Monotone:  Kirchhoff=", monoKirchhoff, " | NonLinearResidual=", monoNonLinear];
        Print[];
    ],
    {caseKey, Keys[$TestCases]}
];

Print[];
Print["=== VERIFICATION SUMMARY ==="];
Print[];
Print["Case                | α    | NL Kirchh  | NL NonLin  | Mono Kirchh | Mono NonLin | NL✓ | Mono✓"];
Print["--------------------|------|------------|------------|-------------|-------------|-----|------"];
Do[
    Print[
        StringPadRight[v["Case"], 20], "| ",
        StringPadRight[ToString[v["Alpha"]], 5], "| ",
        StringPadRight[ToString[NumberForm[v["NL_Kirchhoff"], 4]], 11], "| ",
        StringPadRight[ToString[NumberForm[v["NL_NonLinear"], 4]], 11], "| ",
        StringPadRight[ToString[NumberForm[v["Mono_Kirchhoff"], 4]], 12], "| ",
        StringPadRight[ToString[NumberForm[v["Mono_NonLinear"], 4]], 12], "| ",
        StringPadRight[If[v["NL_Solves_NonLinear"], "YES", "NO"], 4], "| ",
        If[v["Mono_Solves_NonLinear"], "YES", "NO"]
    ],
    {v, allVerifications}
];

Print[];
Print["LEGEND:"];
Print["  Kirchh    = Kirchhoff (flow balance) residual: |KM·j - B|_∞"];
Print["  NonLin    = Non-linear MFG residual: |Nlhs - Nrhs|_∞"];
Print["  NL✓/Mono✓ = Solver's solution actually satisfies the non-linear system (residual < 10⁻⁴)"];
Print[];

(* Export results *)
timestamp = StringReplace[DateString["ISODateTime"], ":" -> "-"];
resultsDir = FileNameJoin[{Directory[], "Results"}];
csvPath = FileNameJoin[{resultsDir, "verify_nonlinear_" <> timestamp <> ".csv"}];
Export[csvPath, allVerifications, "CSV"];
Print["Results exported to: ", csvPath];

Print["Verification complete."];
