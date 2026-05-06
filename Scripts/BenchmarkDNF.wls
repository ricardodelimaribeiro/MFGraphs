#!/usr/bin/env wolframscript

(* Scripts/BenchmarkDNF.wls *)
(* Benchmark comparison between recursive dnfReduce and stack-based dnfReduceProcedural. *)

Needs["MFGraphs`"];

Print["--- DNF Reduction Benchmark ---"];

(* Helper to generate a deep Or-chain of the form:
   And[x1 == 1 || x1 == 0, x2 == 1 || x2 == 0, ..., xn == 1 || xn == 0] *)
generateDeepOrChain[n_] := And @@ Table[Or[Symbol["x" <> ToString[k]] == 0, Symbol["x" <> ToString[k]] == 1], {k, 1, n}];

vars[n_] := Table[Symbol["x" <> ToString[k]], {k, 1, n}];

benchmark[n_] := Module[{expr, allVars, tRec, tProc, resRec, resProc},
    expr = generateDeepOrChain[n];
    allVars = vars[n];
    
    Print["\nTesting n = ", n, " (2^", n, " = ", 2^n, " branches)"];
    
    tRec = AbsoluteTiming[
        resRec = Quiet @ TimeConstrained[solversTools`dnfReduce[True, expr], 10, $Aborted];
    ];
    
    If[MatchQ[resRec, $Aborted | _TerminatedEvaluation],
        Print["  Recursive: Aborted (Timeout or Limit)"],
        Print["  Recursive: ", First[tRec], "s"]
    ];

    tProc = AbsoluteTiming[
        resProc = Quiet @ TimeConstrained[solversTools`Private`dnfReduceProcedural[True, expr, allVars], 10, $Aborted];
    ];
    
    If[MatchQ[resProc, $Aborted | _TerminatedEvaluation],
        Print["  Procedural: Aborted (Timeout or Limit)"],
        Print["  Procedural: ", First[tProc], "s"]
    ];

    If[resRec =!= $Aborted && resProc =!= $Aborted,
        Print["  Results match: ", resRec === resProc]
    ];
];

(* Standard test cases *)
benchmark[4];
benchmark[8];
benchmark[10];

(* Case that may hit $RecursionLimit on some systems if implemented purely recursively without tail-call optimization *)
Block[{$RecursionLimit = 50},
    Print["\nTesting with $RecursionLimit = 50"];
    benchmark[12];
];

Print["\n--- Benchmark Complete ---"];
