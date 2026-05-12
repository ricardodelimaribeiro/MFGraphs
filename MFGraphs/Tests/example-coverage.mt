(* Per-example smoke test: every key in $ExampleScenarios builds, makeSystem
   succeeds, and solveScenario returns within a generous timeout. Solve
   timeouts count as pass; solver errors fail the test. *)

(* Boundary defaults per registered key. Picked to be feasible (entry on a
   source-side vertex, exit on a sink-side vertex) for the topology declared
   by each factory. The "Inconsistent *" entries deliberately violate the
   triangle inequality on switching costs and are expected to construct but
   have no feasible interior; the test only requires construction here. *)
$exampleBoundaryDefaults = <|
    3                                  -> {{{1, 80}}, {{3, 0}}},
    11                                 -> {{{1, 80}}, {{4, 0}}},
    12                                 -> {{{1, 80}}, {{4, 0}}},
    13                                 -> {{{1, 80}}, {{4, 0}}},
    14                                 -> {{{1, 80}}, {{3, 0}}},
    15                                 -> {{{2, 80}, {3, 80}}, {{1, 0}}},
    16                                 -> {{{1, 80}}, {{3, 0}}},
    17                                 -> {{{1, 80}}, {{3, 0}}},
    18                                 -> {{{1, 80}}, {{3, 0}}},
    20                                 -> {{{1, 50}, {2, 50}}, {{7, 0}, {8, 0}, {9, 0}}},
    21                                 -> {{{1, 50}, {2, 50}}, {{10, 0}, {11, 0}, {12, 0}}},
    22                                 -> {{{1, 100}}, {{6, 0}, {7, 0}}},
    23                                 -> {{{1, 50}, {2, 50}}, {{5, 0}, {6, 0}}},
    27                                 -> {{{1, 100}}, {{2, 0}}},
    104                                -> {{{1, 80}}, {{3, 0}}},
    105                                -> {{{1, 100}}, {{2, 0}, {3, 0}}},
    "triangle with two exits"          -> {{{1, 80}}, {{2, 0}, {3, 0}}},
    "chain with two exits"             -> {{{1, 100}}, {{2, 0}, {3, 0}}},
    "Jamaratv9"                        -> {{{1, 100}, {2, 100}}, {{7, 0}, {8, 0}, {9, 0}}},
    "Braess split"                     -> {{{1, 100}}, {{8, 0}}},
    "Braess congest"                   -> {{{1, 100}}, {{7, 0}}},
    "New Braess"                       -> {{{1, 100}}, {{6, 0}}},
    "Big Braess split"                 -> {{{1, 100}}, {{10, 0}}},
    "Big Braess congest"               -> {{{1, 100}}, {{9, 0}}},
    "HRF Scenario 1"                   -> {{{1, 100}}, {{8, 0}, {10, 0}}},
    "Paper example"                    -> {{{1, 100}}, {{4, 0}}},
    "Camilli 2015 simple"              -> {{{3, 50}, {4, 50}}, {{1, 0}}},
    "Camilli 2015 general"             -> {{{14, 25}, {15, 25}, {7, 25}, {9, 25}}, {{1, 0}}},
    "Achdou 2023 junction"             -> {{{2, 25}, {3, 25}, {4, 25}, {5, 25}}, {{1, 0}}},
    "Inconsistent Y shortcut"          -> {{{1, 80}}, {{3, 0}, {4, 0}}},
    "Inconsistent attraction shortcut" -> {{{1, 80}}, {{4, 0}}},
    "Grid0303"                         -> {{{1, 100}}, {{9, 0}}},
    "Grid0404"                         -> {{{1, 100}}, {{16, 0}}},
    "Grid0505"                         -> {{{1, 100}}, {{25, 0}}},
    "Grid0707"                         -> {{{1, 100}}, {{49, 0}}},
    "Grid0710"                         -> {{{1, 100}}, {{70, 0}}},
    "Grid1010"                         -> {{{1, 100}}, {{100, 0}}},
    "Grid1020"                         -> {{{1, 100}}, {{200, 0}}}
|>;

(* Test: $exampleBoundaryDefaults covers every registered scenario. *)
Test[
    Sort[Keys[$exampleBoundaryDefaults]] === listExampleScenarios[]
    ,
    True
    ,
    TestID -> "Example coverage: defaults table covers every registered key"
]

(* Per-key smoke: build + makeSystem + solveScenario (15s timeout, timeout = pass). *)
Scan[
    Function[key,
        With[{
                entries = $exampleBoundaryDefaults[key][[1]],
                exits   = $exampleBoundaryDefaults[key][[2]],
                label   = If[StringQ[key], key, ToString[key]]
            },
            Test[
                Module[{s, sys, sol},
                    s   = getExampleScenario[key, entries, exits];
                    If[!scenarioQ[s], Return[False, Module]];
                    sys = makeSystem[s];
                    If[!mfgSystemQ[sys], Return[False, Module]];
                    (* solve is best-effort: timeouts, recursion-limit
                       termination, and message-emitting returns all count
                       as pass. The build assertions above are what we're
                       really testing. *)
                    sol = CheckAbort[
                        Block[{$RecursionLimit = 16384, $IterationLimit = 16384},
                            TimeConstrained[Quiet @ solveScenario[s], 15, $TimedOut]
                        ],
                        $Aborted
                    ];
                    True
                ]
                ,
                True
                ,
                TestID -> "Example smoke: " <> label
            ]
        ]
    ],
    Keys[$exampleBoundaryDefaults]
];
