#!/usr/bin/env wolframscript
(* PR 3 / Phase B of the disjunct-influence study: benchmark the four
   "DisjunctOrdering" choices on activeSetReduceSystem across all 38
   registered scenarios.

   METHODOLOGY (2026-07, post issue #222 fix):
   - Results collected before the #222 fix were invalid — the option was
     inert and all four orderings executed identical code; Lexicographic-
     first run order turned kernel warm-up into apparent Block-* speedups.
   - Per scenario, one UNTIMED warm-up solve (default path) runs before any
     timed run; each (scenario, mode) pair is then timed $reps times with
     the repetitions INTERLEAVED across modes (mode order rotated per
     repetition and per scenario) so slow kernel-state drift cannot favour
     any one mode. The rep and the run position within the scenario are
     recorded in the CSV so residual order effects stay visible. A mode
     that times out is not re-run in later reps of that scenario (each
     retry would burn the full timeout for no information); skipped reps
     are logged to stdout but not written to the CSV.
   - Since the fix, Block-* orderings force the branch-state path while the
     Lexicographic default routes through dnfReduce, so a plain comparison
     confounds path with ordering. The sweep therefore times five modes:
       Lexicographic     — default path (what a user gets by default)
       Lex-BranchState   — default ordering forced onto the branch-state
                           path (private threshold override): the fair
                           within-path baseline for the Block-* modes
       Block-Vertex / Block-Edge / Block-SCC
   - Each OK result is validity-checked (isValidSystemSolution, 10 s cap).

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
$reps = 3;   (* interleaved timed repetitions per (scenario, mode) *)
$modes = {"Lexicographic", "Lex-BranchState", "Block-Vertex", "Block-Edge", "Block-SCC"};

sanitizeKey[k_Integer] := "case_" <> ToString[k];
sanitizeKey[k_String]  := StringReplace[k, " " -> "_"];

modePath[mode_] := If[mode === "Lexicographic", "dnfReduce", "branchState"];

solveWithMode[sys_, mode_] :=
    Switch[mode,
        "Lex-BranchState",
            (* Default ordering, but forced onto the branch-state path so the
               Block-* modes have a same-path baseline. *)
            Block[{solversTools`Private`$optimizedDNFVarThreshold = Infinity},
                activeSetReduceSystem[sys]],
        _,
            activeSetReduceSystem[sys, "DisjunctOrdering" -> mode]
    ];

(* Time one (scenario, mode, rep) run. position = 1-based slot in this
   scenario's full interleaved run sequence (recorded so order effects
   stay visible). *)
runOne[label_, sys_, mode_, rep_, position_] :=
    Module[{t0, res, status, dt, valid},
        t0 = AbsoluteTime[];
        res = Block[{$RecursionLimit = 16384, $IterationLimit = 16384},
            CheckAbort[
                TimeConstrained[solveWithMode[sys, mode], $timeoutSec, $TimedOut],
                $Aborted]];
        dt = AbsoluteTime[] - t0;
        (* Fully-determined solutions come back as a bare rule LIST;
           parametric/residual ones as an Association. Both are successes.
           (The pre-#222 sweep classified list results as "Other" with
           ruleCount 0 — a second defect in the retracted numbers.) *)
        status = Which[
            res === $TimedOut,                 "Timeout",
            res === $Aborted,                  "Aborted",
            FailureQ[res],                     "Failure",
            AssociationQ[res] || ListQ[res],   "OK",
            True,                              "Other"];
        valid = If[status === "OK",
            TimeConstrained[TrueQ[isValidSystemSolution[sys, res]], 10, "checkTimeout"],
            ""];
        <|
            "scenario"  -> label,
            "ordering"  -> mode,
            "status"    -> status,
            "wallSec"   -> N[dt],
            "ruleCount" -> Which[
                AssociationQ[res], Length[Lookup[res, "Rules", {}]],
                ListQ[res],        Length[res],
                True,              0],
            "rep"       -> rep,
            "position"  -> position,
            "path"      -> modePath[mode],
            "valid"     -> valid
        |>
    ];

$rows = {};
$headers = {"scenario", "ordering", "status", "wallSec", "ruleCount",
            "rep", "position", "path", "valid"};
$skipKeys = {"Grid1020"};   (* kernel-crash on this scenario; see PR 2 notes *)

Print["Sweep: ", Length[$ExampleBoundaryDefaults] - Length[$skipKeys],
      " scenarios x ", Length[$modes], " modes x ", $reps,
      " interleaved reps (rotated order, untimed warm-up), timeout ",
      $timeoutSec, "s"];

(* Write the header and append per-row so a kernel crash mid-sweep
   doesn't lose the prior data. *)
WriteString[$outFile, StringRiffle[$headers, ","] <> "\n"];

writeRow[r_] := WriteString[$outFile,
    StringRiffle[ToString /@ (Lookup[r, #, ""] & /@ $headers), ","] <> "\n"];

$scenarioIndex = 0;
Scan[
    Function[k,
        If[MemberQ[$skipKeys, ToString[k]],
            Print["[skip] ", k, " (kernel-crash list)"],
            Module[{ee = $ExampleBoundaryDefaults[k], scenario, sys, label, runModes},
                label = sanitizeKey[k];
                scenario = getExampleScenario[k, ee[[1]], ee[[2]]];
                sys = makeSystem[scenario];
                (* Untimed warm-up: absorb kernel/internal-cache cold start
                   before any timed run. Capped so a hard scenario cannot
                   stall the sweep twice. *)
                Block[{$RecursionLimit = 16384, $IterationLimit = 16384},
                    Quiet @ CheckAbort[
                        TimeConstrained[activeSetReduceSystem[sys], $timeoutSec, $TimedOut],
                        $Aborted]];
                (* Interleave the reps: within each rep, run every mode once,
                   with the mode order rotated per rep AND per scenario so no
                   mode systematically runs first or right after another. *)
                Module[{position = 0, timedOut = {}},
                    Do[
                        runModes = RotateLeft[$modes,
                            Mod[$scenarioIndex + rep - 1, Length[$modes]]];
                        Scan[
                            Function[ord,
                                position++;
                                If[MemberQ[timedOut, ord],
                                    Print[StringPadRight[label, 22], " ",
                                          StringPadRight[ord, 16],
                                          " rep=", rep,
                                          " SKIPPED (timed out in earlier rep)"],
                                    Module[{r = runOne[label, sys, ord, rep, position]},
                                        If[r["status"] === "Timeout",
                                            AppendTo[timedOut, ord]];
                                        AppendTo[$rows, r];
                                        writeRow[r];
                                        Print[StringPadRight[label, 22], " ",
                                              StringPadRight[ord, 16],
                                              " rep=", rep,
                                              " pos=", position,
                                              " status=", r["status"],
                                              " valid=", r["valid"],
                                              " wall=", r["wallSec"]]
                                    ]
                                ]
                            ],
                            runModes
                        ],
                        {rep, $reps}
                    ]
                ];
                $scenarioIndex++
            ]
        ]
    ],
    Keys[$ExampleBoundaryDefaults]
];

Print["Wrote ", $outFile, " (", Length[$rows], " rows)"];
