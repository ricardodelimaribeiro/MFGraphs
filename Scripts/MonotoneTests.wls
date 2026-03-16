#!/usr/bin/env wolframscript
(* MonotoneTests.wls -- Smoke tests for MonotoneSolver *)

Print["MonotoneSolver Smoke Tests"];
Print["==========================="];
Print[];

PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

$pass = 0;
$fail = 0;

check[label_String, test_] := Module[{result},
  result = Quiet @ Check[test, $Failed];
  If[result === True,
    Print["  PASS: ", label]; $pass++,
    Print["  FAIL: ", label, " (got: ", result, ")"]; $fail++
  ];
];

(* --- Test 1: Degenerate case (Case 1, single vertex) --- *)
Print["--- Degenerate case (Case 1) ---"];
Module[{data, d2e, result},
  data = GetExampleData[1] /. {I1 -> 100, U1 -> 0};
  d2e = DataToEquations[data];
  result = MonotoneSolver[d2e];
  check["Degenerate returns Association", AssociationQ[result]];
];

(* --- Test 2: Linear chain (Case 3), result is an Association with non-negative flows --- *)
Print["--- Linear chain (Case 3) ---"];
Block[{V = Function[{x, edge}, 0]},
  Module[{data, d2e, result, vals},
    data = GetExampleData[3] /. {I1 -> 80, U1 -> 0};
    d2e = DataToEquations[data];
    result = MonotoneSolver[d2e];
    check["Result is Association", AssociationQ[result]];
    If[AssociationQ[result],
      vals = Values[result];
      check["All flows numeric", And @@ (NumericQ /@ vals)];
      check["All flows non-negative", And @@ (# >= -10^-6 & /@ vals)];
    ];
  ];
];

(* --- Test 3: Cached vs uncached equivalence (Case 3) --- *)
Print["--- Cached vs uncached equivalence ---"];
Block[{V = Function[{x, edge}, 0]},
  Module[{data, d2e, cached, uncached, S, halfPairs, signedFlows, rules, substituted,
          jj, cachedEdge, uncachedEdge},
    data = GetExampleData[3] /. {I1 -> 80, U1 -> 0};
    d2e = DataToEquations[data];
    cached = MonotoneSolver[d2e, "UseCachedGradient" -> True];
    uncached = MonotoneSolver[d2e, "UseCachedGradient" -> False];
    check["Both return Associations",
      AssociationQ[cached] && AssociationQ[uncached]];
    If[AssociationQ[cached] && AssociationQ[uncached],
      check["Signed edge flows match within tolerance",
        Norm[Values[cached] - Values[uncached], Infinity] < 10^-4];
    ];
  ];
];

(* --- Summary --- *)
Print[];
Print["==========================="];
Print["Passed: ", $pass, "  Failed: ", $fail];
If[$fail > 0, Exit[1]];
