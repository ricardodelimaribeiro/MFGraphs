(* Wolfram Language Test file *)
(* Tests that the lifted edge-cost field in MonotoneSolver is monotone *)
(* This file should FAIL with a zero-field placeholder and PASS with the *)
(* proper Hessian Riemannian Flow implementation. *)

(* Test: Strict monotonicity of the lifted cost field on an asymmetric Y-network *)
Test[
    Block[{V = Function[{x, edge}, 0],
           alpha = Function[edge, 1],
           g = Function[{m, edge}, -1/m^2]},
        Module[{data, d2e, jj, fieldFn, B, K, jVar, x1Rules, x2Rules, x1, x2, f1, f2},
            data = GetExampleData[7] /. {I1 -> 80, U1 -> 0, U2 -> 10};
            d2e = DataToEquations[data];
            {jj, fieldFn} = MFGraphs`Private`BuildMonotoneField[d2e];
            {B, K} = GetKirchhoffMatrix[d2e][[{1, 2}]];
            jVar = Select[jj, MatchQ[#, _[1, 2, 3]] &][[1]];
            x1Rules = First @ FindInstance[
                K . jj == B && And @@ Thread[jj > 0] && jVar == 20, jj, Reals];
            x2Rules = First @ FindInstance[
                K . jj == B && And @@ Thread[jj > 0] && jVar == 60, jj, Reals];
            x1 = N[jj /. x1Rules];
            x2 = N[jj /. x2Rules];
            f1 = fieldFn[x1];
            f2 = fieldFn[x2];
            Chop[(f1 - f2) . (x1 - x2), 10^-8] > 10^-6
        ]
    ]
    ,
    True
    ,
    TestID -> "Monotonicity: (f(x1)-f(x2)).(x1-x2) > 0 on Y-network"
]

(* Test: Field is non-trivial (guard against zero placeholder) *)
Test[
    Block[{V = Function[{x, edge}, 0],
           alpha = Function[edge, 1],
           g = Function[{m, edge}, -1/m^2]},
        Module[{data, d2e, jj, fieldFn, B, K, jVar, x1Rules, x2Rules, x1, x2, f1, f2},
            data = GetExampleData[7] /. {I1 -> 80, U1 -> 0, U2 -> 10};
            d2e = DataToEquations[data];
            {jj, fieldFn} = MFGraphs`Private`BuildMonotoneField[d2e];
            {B, K} = GetKirchhoffMatrix[d2e][[{1, 2}]];
            jVar = Select[jj, MatchQ[#, _[1, 2, 3]] &][[1]];
            x1Rules = First @ FindInstance[
                K . jj == B && And @@ Thread[jj > 0] && jVar == 20, jj, Reals];
            x2Rules = First @ FindInstance[
                K . jj == B && And @@ Thread[jj > 0] && jVar == 60, jj, Reals];
            x1 = N[jj /. x1Rules];
            x2 = N[jj /. x2Rules];
            f1 = fieldFn[x1];
            f2 = fieldFn[x2];
            Norm[f1 - f2, Infinity] > 10^-6
        ]
    ]
    ,
    True
    ,
    TestID -> "Monotonicity: field is non-trivial (not zero placeholder)"
]
