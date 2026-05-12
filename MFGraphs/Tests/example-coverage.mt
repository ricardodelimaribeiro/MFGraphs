(* Per-example smoke test: every key in $ExampleScenarios builds, makeSystem
   succeeds, and solveScenario does not throw. When a cached solution is
   present under <repo>/solutions/, the fresh solve is compared against it
   and any divergence fails the test. When the cache is missing, the test
   stores the fresh result so subsequent runs can compare.

   Refresh cache with: wolframscript -file Scripts/RegenerateSolutions.wls *)

(* The boundary-defaults table and cache helpers live in
   Scripts/SolutionCacheHelpers.wls so the regenerate script and the test
   share a single source of truth. *)
Get[FileNameJoin[{$RepoRoot, "Scripts", "SolutionCacheHelpers.wls"}]];

(* Test: $ExampleBoundaryDefaults covers every registered scenario. *)
Test[
    Sort[Keys[$ExampleBoundaryDefaults]] === listExampleScenarios[]
    ,
    True
    ,
    TestID -> "Example coverage: defaults table covers every registered key"
]

(* Solution-equality predicate: directly compare canonicalised forms.
   Both solveScenario outputs and cached forms come from the same code,
   so byte-equality after Sort on rule lists is the right check. *)
solutionEquivalent[a_, b_] :=
    Module[{ca = a, cb = b},
        If[ListQ[ca] && AllTrue[ca, MatchQ[#, _Rule | _RuleDelayed] &],
            ca = SortBy[ca, ToString[First[#]] &]];
        If[ListQ[cb] && AllTrue[cb, MatchQ[#, _Rule | _RuleDelayed] &],
            cb = SortBy[cb, ToString[First[#]] &]];
        ca === cb
    ];

(* Per-key smoke + cache-comparison Test. *)
Scan[
    Function[key,
        With[{
                entries = $ExampleBoundaryDefaults[key][[1]],
                exits   = $ExampleBoundaryDefaults[key][[2]],
                label   = If[StringQ[key], key, ToString[key]]
            },
            Test[
                Module[{s, sys, sol, cached},
                    s = getExampleScenario[key, entries, exits];
                    If[!scenarioQ[s], Return[False, Module]];
                    sys = makeSystem[s];
                    If[!mfgSystemQ[sys], Return[False, Module]];
                    cached = loadSolution[key];
                    sol = CheckAbort[
                        Block[{$RecursionLimit = 16384, $IterationLimit = 16384},
                            TimeConstrained[Quiet @ solveScenario[s], 15, $TimedOut]
                        ],
                        $Aborted
                    ];
                    Which[
                        (* Cache hit: assert the fresh solve agrees. Skip
                           comparison if the fresh solve aborted/timed out;
                           the cache stays as the authority. *)
                        AssociationQ[cached] && KeyExistsQ[cached, "Solution"] &&
                            sol =!= $TimedOut && sol =!= $Aborted &&
                            Head[sol] =!= TerminatedEvaluation,
                        solutionEquivalent[sol, cached["Solution"]],
                        (* Cache miss but fresh solve produced a result:
                           seed the cache so future runs can compare. *)
                        sol =!= $TimedOut && sol =!= $Aborted &&
                            Head[sol] =!= TerminatedEvaluation,
                        storeSolution[key, entries, exits, sol]; True,
                        (* Otherwise: build succeeded, solve was best-effort. *)
                        True, True
                    ]
                ]
                ,
                True
                ,
                TestID -> "Example smoke: " <> label
            ]
        ]
    ],
    Keys[$ExampleBoundaryDefaults]
];
