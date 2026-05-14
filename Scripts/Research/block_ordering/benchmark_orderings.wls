#!/usr/bin/env wolframscript
(* PR 3 / Phase B of the disjunct-influence study: benchmark the four
   "DisjunctOrdering" choices on activeSetReduceSystem across all 38
   registered scenarios.

   Per (scenario, ordering): wall time, Status, MaxRecursionDepth,
   BranchesStarted/Kept/Dropped (from dnfReduceDiagnosticReport's Summary
   when available; otherwise N/A — activeSetReduceSystem itself does not
   emit those counters, so we time the solve and separately call the
   diagnostic for branching shape).

   Writes one CSV row per (scenario, ordering) to:
     Scripts/Research/results/block_ordering/benchmark.csv

   The benchmark wraps each run in TimeConstrained[..., 60] so a single
   slow ordering can't stall the sweep. *)

$RepoRoot = ParentDirectory[ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]]];
PrependTo[$Path, FileNameJoin[{$RepoRoot, "MFGraphs"}]];
Needs["MFGraphs`"];
Get[FileNameJoin[{$RepoRoot, "Scripts", "SolutionCacheHelpers.wls"}]];

$resultsDir = FileNameJoin[{$RepoRoot, "Scripts", "Research", "results", "block_ordering"}];
If[!DirectoryQ[$resultsDir],
    CreateDirectory[$resultsDir, CreateIntermediateDirectories -> True]];
$outFile = FileNameJoin[{$resultsDir, "benchmark.csv"}];

$timeoutSec = 60;
$orderings = {"Lexicographic", "Block-Vertex", "Block-Edge", "Block-SCC"};

sanitizeKey[k_Integer] := "case_" <> ToString[k];
sanitizeKey[k_String]  := StringReplace[k, " " -> "_"];

(* Time one (scenario, ordering) pair. *)
runOne[label_, sys_, ordering_] :=
    Module[{t0, res, status, dt},
        t0 = AbsoluteTime[];
        res = Block[{$RecursionLimit = 16384, $IterationLimit = 16384},
            CheckAbort[
                TimeConstrained[
                    activeSetReduceSystem[sys, "DisjunctOrdering" -> ordering],
                    $timeoutSec, $TimedOut],
                $Aborted]];
        dt = AbsoluteTime[] - t0;
        status = Which[
            res === $TimedOut, "Timeout",
            res === $Aborted,  "Aborted",
            AssociationQ[res], "OK",
            True,              "Other"];
        <|
            "scenario"  -> label,
            "ordering"  -> ordering,
            "status"    -> status,
            "wallSec"   -> N[dt],
            "ruleCount" -> If[AssociationQ[res], Length[Lookup[res, "Rules", {}]], 0]
        |>
    ];

$rows = {};
$headers = {"scenario", "ordering", "status", "wallSec", "ruleCount"};
$skipKeys = {"Grid1020"};   (* kernel-crash on this scenario; see PR 2 notes *)

Print["Sweep: ", Length[$ExampleBoundaryDefaults] - Length[$skipKeys],
      " scenarios x ", Length[$orderings], " orderings, timeout ", $timeoutSec, "s"];

(* Write the header and append per-row so a kernel crash mid-sweep
   doesn't lose the prior data. *)
WriteString[$outFile, StringRiffle[$headers, ","] <> "\n"];

writeRow[r_] := WriteString[$outFile,
    StringRiffle[ToString /@ (Lookup[r, #, ""] & /@ $headers), ","] <> "\n"];

Scan[
    Function[k,
        If[MemberQ[$skipKeys, ToString[k]],
            Print["[skip] ", k, " (kernel-crash list)"],
            Module[{ee = $ExampleBoundaryDefaults[k], scenario, sys, label},
                label = sanitizeKey[k];
                scenario = getExampleScenario[k, ee[[1]], ee[[2]]];
                sys = makeSystem[scenario];
                Scan[
                    Function[ord,
                        Module[{r = runOne[label, sys, ord]},
                            AppendTo[$rows, r];
                            writeRow[r];
                            Print[StringPadRight[label, 22], " ",
                                  StringPadRight[ord, 14],
                                  " status=", r["status"], " wall=", r["wallSec"]]
                        ]
                    ],
                    $orderings
                ]
            ]
        ]
    ],
    Keys[$ExampleBoundaryDefaults]
];

Print["Wrote ", $outFile, " (", Length[$rows], " rows)"];
