#!/usr/bin/env wolframscript
(* PR 1 pilot: measure per-disjunct branching counts on 3 representative
   scenarios. Validates that the new BranchesStartedByDisjunct /
   BranchesKeptByDisjunct / BranchesDroppedByDisjunct counters introduced
   in PR 1 of the disjunct-influence study work as intended.

   Run:
     wolframscript -file Scripts/Research/disjunct_influence/measure_branching.wls

   Writes one JSON per scenario to:
     Scripts/Research/results/branching/<scenario>.json

   Each JSON has shape:
     {
       "scenario": "<key>",
       "totals": {"BranchesStarted": <int>, "BranchesKept": <int>, "BranchesDropped": <int>},
       "perDisjunct": [
         {"disjunctId": "<canonical>", "started": <int>, "kept": <int>, "dropped": <int>,
          "pruneRate": <real>, "survivalRate": <real>},
         ...
       ]
     }

   Sorted descending by `started` so the most-branching disjuncts appear first. *)

$RepoRoot = ParentDirectory[ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]]];
PrependTo[$Path, FileNameJoin[{$RepoRoot, "MFGraphs"}]];
Needs["MFGraphs`"];
Get[FileNameJoin[{$RepoRoot, "Scripts", "SolutionCacheHelpers.wls"}]];

$resultsDir = FileNameJoin[{$RepoRoot, "Scripts", "Research", "results", "branching"}];
If[!DirectoryQ[$resultsDir], CreateDirectory[$resultsDir, CreateIntermediateDirectories -> True]];

(* Build the full sweep from $ExampleBoundaryDefaults (the canonical
   per-key entries/exits used by example-coverage.mt and the cache).
   sanitizeKey turns numeric keys (e.g. 22) and string keys (e.g.
   "Jamaratv9") into safe filenames. *)
sanitizeKey[k_Integer] := "case_" <> ToString[k];
sanitizeKey[k_String]  := StringReplace[k, " " -> "_"];

$sweep = Function[k,
    Module[{ee = $ExampleBoundaryDefaults[k]},
        {sanitizeKey[k], k, ee[[1]], ee[[2]]}
    ]
] /@ Keys[$ExampleBoundaryDefaults];

$timeoutSec = 60;

measureScenario[label_String, scenario_] :=
    Module[{sys, t0, rep, sm, started, kept, dropped, perDisj, sorted, json,
            outFile, conservation},
        sys = makeSystem[scenario];
        t0  = AbsoluteTime[];
        rep = Block[{$RecursionLimit = 16384, $IterationLimit = 16384},
            CheckAbort[
                dnfReduceDiagnosticReport[sys, "Timeout" -> $timeoutSec],
                <|"Status" -> "Aborted",
                  "Summary" -> <|"BranchesStarted" -> 0, "BranchesKept" -> 0,
                                 "BranchesDropped" -> 0, "MaxRecursionDepth" -> 0,
                                 "BranchesStartedByDisjunct" -> <||>,
                                 "BranchesKeptByDisjunct" -> <||>,
                                 "BranchesDroppedByDisjunct" -> <||>|>|>
            ]
        ];
        sm  = rep["Summary"];
        started = sm["BranchesStartedByDisjunct"];
        kept    = sm["BranchesKeptByDisjunct"];
        dropped = sm["BranchesDroppedByDisjunct"];
        (* BranchesDropped includes drops from non-Or-branch sources (the
           And-block elem===False case in dnfReduceInstrumentedAnd) which
           have no associated disjunct. The gap (BranchesDropped -
           sum(DroppedByDisjunct)) is the count of those non-Or drops.
           BranchesStarted and BranchesKept always come from Or contexts so
           they should conserve exactly. *)
        conservation = <|
            "BranchesStarted" -> sm["BranchesStarted"] === Total[Values[started]],
            "BranchesKept"    -> sm["BranchesKept"]    === Total[Values[kept]],
            "DropGapNonOr"    -> sm["BranchesDropped"] - Total[Values[dropped]]
        |>;
        perDisj = Function[id,
            Module[{s = Lookup[started, id, 0], k = Lookup[kept, id, 0],
                    d = Lookup[dropped, id, 0]},
                <|
                    "disjunctId"  -> id,
                    "started"     -> s,
                    "kept"        -> k,
                    "dropped"     -> d,
                    "pruneRate"   -> N[d / Max[1, s]],
                    "survivalRate" -> N[k / Max[1, s]]
                |>
            ]
        ] /@ Keys[started];
        sorted = Reverse @ SortBy[perDisj, #["started"] &];
        json = <|
            "scenario"     -> label,
            "elapsedSec"   -> N[AbsoluteTime[] - t0],
            "status"       -> rep["Status"],
            "totals"       -> <|
                "BranchesStarted"   -> sm["BranchesStarted"],
                "BranchesKept"      -> sm["BranchesKept"],
                "BranchesDropped"   -> sm["BranchesDropped"],
                "UniqueDisjuncts"   -> Length[Keys[started]],
                "MaxRecursionDepth" -> sm["MaxRecursionDepth"]
            |>,
            "conservation" -> conservation,
            "perDisjunct"  -> sorted
        |>;
        outFile = FileNameJoin[{$resultsDir, label <> ".json"}];
        Export[outFile, json, "RawJSON"];
        Print[
            StringPadRight[label, 12],
            " status=", rep["Status"],
            " started=", sm["BranchesStarted"],
            " unique=", Length[Keys[started]],
            " conservation=", conservation["BranchesStarted"] && conservation["BranchesKept"],
            " dropGap=", conservation["DropGapNonOr"],
            " -> ", FileNameTake[outFile]
        ]
    ];

Print["Sweep: ", Length[$sweep], " scenarios, timeout ", $timeoutSec, "s, results to ", $resultsDir];
Scan[
    Function[entry,
        Module[{label, key, entries, exits, scenario},
            {label, key, entries, exits} = entry;
            scenario = getExampleScenario[key, entries, exits];
            measureScenario[label, scenario]
        ]
    ],
    $sweep
];
Print["Done."];
