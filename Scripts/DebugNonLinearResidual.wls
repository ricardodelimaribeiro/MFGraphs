#!/usr/bin/env wolframscript
(* DebugNonLinearResidual.wls

   Detailed inspection of why NonLinearSolver/MonotoneSolver solutions
   don't satisfy the non-linear MFG equations Nlhs == Nrhs.
*)

SetDirectory[DirectoryName[$ScriptCommandLine[[1]]] <> "/.."];
Get["MFGraphs/MFGraphs.wl"];

Print["NonLinear Residual Debug"];
Print["========================"];

(* Test on Y7 alpha=2 *)
data = GetExampleData[7] /. {I1 -> 50, U1 -> 0, U2 -> 0};
d2e = DataToEquations[data];
alphaFn = Function[edge, 2];

Print["d2e keys: ", Keys[d2e]];
Print[];

Nlhs = Lookup[d2e, "Nlhs"];
Nrhs = Lookup[d2e, "Nrhs"];

Print["Nlhs head: ", Head[Nlhs], ", length: ", If[ListQ[Nlhs], Length[Nlhs], "scalar"]];
Print["Nrhs head: ", Head[Nrhs], ", length: ", If[ListQ[Nrhs], Length[Nrhs], "scalar"]];
Print[];

Print["First 3 Nlhs entries:"];
Do[Print["  [", i, "] ", Nlhs[[i]]], {i, 1, Min[3, Length[Nlhs]]}];
Print[];
Print["First 3 Nrhs entries:"];
Do[Print["  [", i, "] ", Nrhs[[i]]], {i, 1, Min[3, Length[Nrhs]]}];
Print[];

Print["Variables in Nlhs[[1]]: ", Variables[Nlhs[[1]]]];
Print["Variables in Nrhs[[1]]: ", Variables[Nrhs[[1]]]];
Print[];

(* Run NonLinearSolver *)
Print["Running NonLinearSolver..."];
nlResult = Quiet @ NonLinearSolver[d2e,
    "CongestionExponentFunction" -> alphaFn,
    "MaxIterations" -> 15, "Tolerance" -> 10^-6];

nlSol = Lookup[nlResult, "AssoNonCritical"];
Print["NL solution keys (first 10): ", Take[Keys[nlSol], Min[10, Length[nlSol]]]];
Print["NL u variables: ",
  Select[Normal[nlSol], MatchQ[First[#], _MFGraphs`Private`u | u[__]] &]
];
Print[];
Print["NL j variables (first 5):"];
Do[
  Print["  ", n[[1]], " -> ", n[[2]]],
  {n, Take[Select[Normal[nlSol], MatchQ[First[#], _MFGraphs`Private`j] &], Min[5, Length[nlSol]]]}
];
Print[];

(* Now evaluate Nlhs and Nrhs *)
Print["Evaluating Nlhs /. nlSol on first entry:"];
Print["  Nlhs[[1]] = ", Nlhs[[1]]];
Print["  After substitution: ", Nlhs[[1]] /. Normal[nlSol]];
Print["  Numeric value: ", N[Nlhs[[1]] /. Normal[nlSol]]];
Print[];
Print["  Nrhs[[1]] = ", Nrhs[[1]]];
Print["  After substitution: ", Nrhs[[1]] /. Normal[nlSol]];
Print["  Numeric value: ", N[Nrhs[[1]] /. Normal[nlSol]]];
Print[];

(* Show the worst residual *)
lhsVals = N[Nlhs /. Normal[nlSol]];
rhsVals = N[Nrhs /. Normal[nlSol]];
diffs = lhsVals - rhsVals;
Print["All residuals (lhs - rhs):"];
Do[Print["  [", i, "] ", diffs[[i]], "  (lhs=", lhsVals[[i]], ", rhs=", rhsVals[[i]], ")"],
   {i, 1, Length[diffs]}];
Print[];
Print["Max abs residual: ", Max[Abs[diffs]]];
Print[];

(* Check the IsNonLinearSolution validation function *)
Print["IsNonLinearSolution check:"];
isValid = IsNonLinearSolution[nlResult];
Print["  Result: ", isValid];
