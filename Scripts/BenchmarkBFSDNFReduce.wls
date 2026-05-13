#!/usr/bin/env wolframscript
(* Head-to-head: dnfReduceSystem (DFS) vs bfsDNFReduceSystem (BFS) on a few
   reference scenarios plus case 21. Reports elapsed time, validity, and
   cachedSolve hit/miss counts.

   Usage: wolframscript -file Scripts/BenchmarkBFSDNFReduce.wls
                       [--timeout s]   default 60
                       [--case name]   default all  (one of: grid-2x3, grid-3x2, grid-3x3, case-12, case-21) *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength   = 0;

$Timeout    = 60;
$CaseFilter = "all";

Module[{args = Rest[$ScriptCommandLine], i = 1, arg},
  While[i <= Length[args],
    arg = args[[i]];
    Switch[arg,
      "--timeout",
        i++; If[i <= Length[args],
          Module[{v = Quiet @ Check[ToExpression[args[[i]]], $Failed]},
            If[NumericQ[v] && v > 0, $Timeout = v]]],
      "--case",
        i++; If[i <= Length[args], $CaseFilter = ToLowerCase[args[[i]]]],
      _, Null
    ];
    i++
  ]
];

$Cases = {
  {"grid-2x3", Function[gridScenario[{2, 3}, {{1, 100}}, {{6, 0}}]]},
  {"grid-3x2", Function[gridScenario[{3, 2}, {{1, 100}}, {{6, 0}}]]},
  {"grid-3x3", Function[gridScenario[{3, 3}, {{1, 100}}, {{9, 0}}]]},
  {"case-12",  Function[getExampleScenario[12, {{1, 100}}, {{4, 0}}]]},
  {"case-21",  Function[getExampleScenario[21, {{1, 50}, {2, 50}}, {{10, 0}, {11, 0}, {12, 0}}]]}
};

If[$CaseFilter =!= "all",
  $Cases = Select[$Cases, #[[1]] === $CaseFilter &];
  If[$Cases === {}, Print["No case named '", $CaseFilter, "'."]; Exit[2]]
];

Print["MFGraphs DFS-vs-BFS DNF reducer benchmark"];
Print["Date\t",    DateString[]];
Print["Timeout\t", $Timeout, " s"];
Print[];
Print[StringRiffle[
  {"Case", "Variant", "Elapsed_s", "Status", "Valid"},
  "\t"]];

Fmt[x_?NumericQ] := ToString[NumberForm[N[x], {Infinity, 4}, NumberPadding -> {"", "0"}]];
Fmt[x_]          := ToString[x];

runOne[name_String, buildFn_, variant_String] :=
    Module[{sys, t, res, status, valid},
        sys = makeSystem[buildFn[]];
        {t, res} = AbsoluteTiming @ TimeConstrained[
            Switch[variant,
                "DFS", dnfReduceSystem[sys],
                "BFS", bfsDNFReduceSystem[sys]
            ],
            $Timeout, $TimedOut
        ];
        status = Which[
            res === $TimedOut, "Timeout",
            FailureQ[res],     "Failure",
            True,              "OK"
        ];
        valid = If[status === "OK", isValidSystemSolution[sys, res], "n/a"];
        Print[StringRiffle[{name, variant, Fmt[t], status, Fmt[valid]}, "\t"]];
    ];

Scan[
    Function[case,
        runOne[case[[1]], case[[2]], "DFS"];
        runOne[case[[1]], case[[2]], "BFS"];
    ],
    $Cases
];

(* :!CodeAnalysis::EndBlock:: *)
