#!/usr/bin/env wolframscript
(* Sweep $optimizedDNFVarThreshold across multiple values to find the empirical
   crossover between branch-state enumeration and native dnfReduce. Mutates the
   private symbol at runtime — no source edit. Skips dnfReduceSystem (doesn't
   read the threshold). *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

Print["MFGraphs $optimizedDNFVarThreshold sweep"];
Print["======================================="];
Print["Date: ", DateString[]];
Print[];

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength = 0;

$Thresholds = {0, 4, 6, 8, 10, 12, 16, Infinity};
$Solvers = {
  {"optimizeddnf", solversTools`optimizedDNFReduceSystem},
  {"activeset",    solversTools`activeSetReduceSystem}
};

$BenchCases = {
  {"chain-2v",         Function[gridScenario[{2}, {{1, 10}}, {{2, 0}}]]},
  {"chain-3v-1exit",   Function[gridScenario[{3}, {{1, 10}}, {{3, 0}}]]},
  {"chain-3v-2exit",   Function[gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}]]},
  {"chain-3-midentry", Function[gridScenario[{3}, {{2, 100.0}}, {{1, 0.0}, {3, 10.0}}]]},
  {"chain-5v-1exit",   Function[gridScenario[{5}, {{1, 10}}, {{5, 0}}]]},
  {"grid-2x3",         Function[gridScenario[{2, 3}, {{1, 100}}, {{6, 0}}]]},
  {"grid-3x2",         Function[gridScenario[{3, 2}, {{1, 100}}, {{6, 0}}]]},
  {"example-12",       Function[getExampleScenario[12, {{1, 100.0}}, {{4, 0.0}}]]}
};

$TimeoutSecs = 30;

(* allVars count exposed via buildSolverInputs (same path the solvers take). *)
allVarsCount[sys_] :=
  Module[{inputs},
    inputs = solversTools`Private`buildSolverInputs[sys];
    If[ListQ[inputs] && Length[inputs] >= 2, Length[inputs[[2]]], Missing["NotAvailable"]]
  ];

(* RepeatedTiming auto-tunes sample count to a stable mean. *)
repeatMs[fn_] :=
  Module[{t, result},
    {t, result} = RepeatedTiming[TimeConstrained[fn[], $TimeoutSecs, $TimedOut]];
    If[result === $TimedOut, Missing["TimedOut"], Round[N[t * 1000], 0.01]]
  ];

(* Precompute systems once; threshold doesn't affect makeSystem. *)
Print["Building systems..."];
$Systems = AssociationThread[
  $BenchCases[[All, 1]],
  Function[{c}, makeSystem[c[[2]][]]] /@ $BenchCases
];
$AllVars = AssociationThread[
  $BenchCases[[All, 1]],
  allVarsCount /@ Values[$Systems]
];
Print["allVars per case: ", $AllVars];
Print[];

(* Warm up the kernel: each (system, solver) once at default threshold so
   first-iteration JIT cost doesn't bias one cell. *)
Print["Warming up..."];
solversTools`Private`$optimizedDNFVarThreshold = 8;
Do[
  Do[
    TimeConstrained[solver[[2]][$Systems[caseName]], $TimeoutSecs, $TimedOut],
    {solver, $Solvers}
  ],
  {caseName, $BenchCases[[All, 1]]}
];
Print[];

$Results = Reap[
  Do[
    solversTools`Private`$optimizedDNFVarThreshold = thr;
    confirmedThr = solversTools`Private`$optimizedDNFVarThreshold;
    If[confirmedThr =!= thr,
      Print["WARN: threshold mutation did not stick: requested ", thr, " got ", confirmedThr]
    ];
    Do[
      Do[
        caseName  = case[[1]];
        sys       = $Systems[caseName];
        solverFn  = solver[[2]];
        repMs     = repeatMs[Function[Null, solverFn[sys]]];
        Print[
          "thr=", thr, " | ", case[[1]], " | ", solver[[1]],
          " | RepMs=", repMs
        ];
        Sow[<|
          "Threshold" -> thr,
          "Case"      -> caseName,
          "Solver"    -> solver[[1]],
          "AllVars"   -> $AllVars[caseName],
          "RepMs"     -> repMs
        |>];
        ,
        {solver, $Solvers}
      ];
      ,
      {case, $BenchCases}
    ];
    ,
    {thr, $Thresholds}
  ];
][[2, 1]];

Print[];
Print["=== Pivot table (RepeatedTiming average ms) ==="];

(* Build per (case, solver) row across thresholds. *)
$Rows = GroupBy[$Results, {#Case, #Solver} &];
thrLabel[thr_] := If[thr === Infinity, "Inf", ToString[thr]];

PrintRow[label_, cells_List] :=
  Print["| ", StringRiffle[ToString /@ Flatten[{label, cells}], " | "], " |"];

Print["| Case | allVars | Solver | thr=", StringRiffle[thrLabel /@ $Thresholds, " | thr="], " |"];
Print["|------|--------:|--------|", StringJoin[Table["------:|", Length[$Thresholds]]]];

KeyValueMap[
  Function[{key, rows},
    cells = Table[
      With[{r = SelectFirst[rows, #Threshold === thr &, <||>]},
        ToString @ Lookup[r, "RepMs", "-"]
      ],
      {thr, $Thresholds}
    ];
    PrintRow[
      {key[[1]], First[rows]["AllVars"], key[[2]]},
      cells
    ];
  ],
  $Rows
];

(* Per-(case, solver) optimum — minimal RepMs threshold value. *)
Print[];
Print["=== Per (case, solver) optimum threshold ==="];
optimumPerRow = KeyValueMap[
  Function[{key, rows},
    With[{numeric = Cases[rows, r_ /; NumericQ[r["RepMs"]] :> r]},
      If[numeric === {},
        {key[[1]], key[[2]], "—", "—"},
        With[{best = MinimalBy[numeric, #RepMs &][[1]]},
          {key[[1]], key[[2]], thrLabel[best["Threshold"]], best["RepMs"]}
        ]
      ]
    ]
  ],
  $Rows
];
Print["| Case | Solver | OptThreshold | OptRepMs |"];
Print["|------|--------|--------------|---------:|"];
(PrintRow[#[[1]], Rest[#]] &) /@ optimumPerRow;

(* Safe global threshold: max of per-row optimum thresholds (excluding rows with
   no numeric optimum and excluding sub-millisecond rows that sit in noise). *)
numericOpts = Cases[
  optimumPerRow,
  {_, _, t_String, ms_?NumericQ} /; ms >= 1.0 :> If[t === "Inf", Infinity, ToExpression[t]]
];
safeGlobal = If[numericOpts === {}, Missing["NoData"], Max[numericOpts]];
Print[];
Print["Safe global threshold (max-of-per-row, excluding sub-ms noise): ", safeGlobal];
Print["Default value: 0 (set in solversTools.wl). Sweep values: ", $Thresholds];

$ResultsDir = FileNameJoin[{$RepoRoot, "Results"}];
If[!DirectoryQ[$ResultsDir], CreateDirectory[$ResultsDir]];
$Timestamp = DateString[{"Year", "Month", "Day", "-", "Hour", "Minute", "Second"}];
$CSVFile = FileNameJoin[{$ResultsDir, "branch_state_threshold_sweep_" <> $Timestamp <> ".csv"}];
Export[
  $CSVFile,
  Prepend[
    Map[{#Threshold /. Infinity -> "Infinity", #Case, #Solver, #AllVars, #RepMs} &, $Results],
    {"Threshold", "Case", "Solver", "AllVars", "RepMs"}
  ],
  "CSV"
];
Print[];
Print["CSV written to: ", $CSVFile];

(* :!CodeAnalysis::EndBlock:: *)
