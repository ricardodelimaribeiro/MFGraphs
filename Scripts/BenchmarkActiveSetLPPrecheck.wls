#!/usr/bin/env wolframscript
(* Head-to-head: activeSetReduceSystem with and without LPPrecheck on a few
   reference scenarios plus case 21. Reports elapsed time, validity,
   precheck-call count, and prune count.

   Usage: wolframscript -file Scripts/BenchmarkActiveSetLPPrecheck.wls
                       [--timeout s]      default 60
                       [--lp-timeout s]   default 0.2
                       [--case name]      default all  (one of: grid-2x3, grid-3x2, grid-3x3, case-12, case-21) *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength = 0;

$Timeout    = 60;
$LpTimeout  = 0.2;
$CaseFilter = "all";

Module[{args = Rest[$ScriptCommandLine], i = 1, arg},
  While[i <= Length[args],
    arg = args[[i]];
    Switch[arg,
      "--timeout",
        i++; If[i <= Length[args],
          Module[{v = Quiet @ Check[ToExpression[args[[i]]], $Failed]},
            If[NumericQ[v] && v > 0, $Timeout = v]]],
      "--lp-timeout",
        i++; If[i <= Length[args],
          Module[{v = Quiet @ Check[ToExpression[args[[i]]], $Failed]},
            If[NumericQ[v] && v > 0, $LpTimeout = v]]],
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
  If[$Cases === {},
    Print["No case named '", $CaseFilter, "'."];
    Exit[2]
  ]
];

Print["MFGraphs activeSetReduceSystem LP-precheck benchmark"];
Print["Date\t", DateString[]];
Print["Timeout\t", $Timeout, " s"];
Print["LP-timeout\t", $LpTimeout, " s/call"];
Print[];
Print[StringRiffle[
  {"Case", "Variant", "Elapsed_s", "Status", "Valid", "PrecheckCalls", "Pruned", "PruneRate"},
  "\t"]];

Fmt[x_?NumericQ] := ToString[NumberForm[N[x], {Infinity, 4}, NumberPadding -> {"", "0"}]];
Fmt[x_]          := ToString[x];

runOne[name_String, buildFn_, useLP_] :=
  Module[{sys, t, res, calls, pruned, status, valid},
    sys = makeSystem[buildFn[]];
    {t, res} = AbsoluteTiming @ TimeConstrained[
      If[useLP,
        activeSetReduceSystem[sys, "LPPrecheck" -> True, "LPPrecheckTimeout" -> $LpTimeout],
        activeSetReduceSystem[sys]
      ],
      $Timeout, $TimedOut
    ];
    calls  = solversTools`Private`$lpPrecheckCalls;
    pruned = solversTools`Private`$lpPrecheckPruned;
    status = Which[
      res === $TimedOut,         "Timeout",
      FailureQ[res],             "Failure",
      True,                      "OK"
    ];
    valid = If[status === "OK", isValidSystemSolution[sys, res], "n/a"];
    Print[StringRiffle[
      {name,
       If[useLP, "LPPrecheck=on", "LPPrecheck=off"],
       Fmt[t],
       status,
       Fmt[valid],
       Fmt[calls],
       Fmt[pruned],
       Fmt[N[pruned/Max[1, calls]]]},
      "\t"]];
  ];

Scan[
  Function[case,
    runOne[case[[1]], case[[2]], False];
    runOne[case[[1]], case[[2]], True];
  ],
  $Cases
];

(* :!CodeAnalysis::EndBlock:: *)
