#!/usr/bin/env wolframscript
(* Regenerate cached solutions for the registered example scenarios.
   Solutions are written under <repo>/solutions/<sanitized-key>.wxf.

   Usage:
     wolframscript -file Scripts/RegenerateSolutions.wls
     wolframscript -file Scripts/RegenerateSolutions.wls Jamaratv9 Grid0303
     wolframscript -file Scripts/RegenerateSolutions.wls --timeout 600

   Scenarios that exceed the per-scenario timeout are skipped (no file
   written; existing cache file, if any, is left in place). *)

$RepoRoot = ParentDirectory[DirectoryName[$InputFileName]];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;  (* overridden below if --verbose set *)

Get[FileNameJoin[{$RepoRoot, "Scripts", "SolutionCacheHelpers.wls"}]];

(* --- argument parsing --- *)

$rawArgs = $ScriptCommandLine[[2 ;;]];

$timeoutSec = 300;
$selectedKeys = {};
$verbose = False;

processArgs[args_List] :=
    Module[{rest = args, head},
        While[rest =!= {},
            head = First[rest];
            Which[
                head === "--timeout" || head === "-t",
                If[Length[rest] < 2,
                    Print["Missing value for --timeout"]; Exit[2]];
                $timeoutSec = ToExpression[rest[[2]]];
                rest = rest[[3 ;;]],
                head === "--verbose" || head === "-v",
                $verbose = True;
                rest = Rest[rest],
                StringStartsQ[head, "--"],
                Print["Unknown flag: ", head]; Exit[2],
                True,
                AppendTo[$selectedKeys, head];
                rest = Rest[rest]
            ]
        ];
    ];

processArgs[$rawArgs];

If[$verbose, $MFGraphsVerbose = True];

If[$selectedKeys === {}, $selectedKeys = Keys[$ExampleBoundaryDefaults]];

(* String args may match either string keys or numeric keys (case_3 → 3). *)
resolveKey[arg_] :=
    Module[{numericMatch, stringMatch},
        If[KeyExistsQ[$ExampleBoundaryDefaults, arg],
            Return[arg]];
        numericMatch = Quiet @ ToExpression[arg];
        If[IntegerQ[numericMatch] && KeyExistsQ[$ExampleBoundaryDefaults, numericMatch],
            Return[numericMatch]];
        stringMatch = arg;
        If[KeyExistsQ[$ExampleBoundaryDefaults, stringMatch],
            Return[stringMatch]];
        Print["Unknown scenario key: ", arg];
        Exit[3]
    ];

$selectedKeys = If[Head[#] === String, resolveKey[#], #] & /@ $selectedKeys;

Print["Regenerating ", Length[$selectedKeys], " scenario(s) with timeout ", $timeoutSec, "s."];
Print["Cache dir: ", solutionCacheDir[]];

(* --- main loop --- *)

results = {};

Scan[
    Function[key,
        Module[{ee, entries, exits, s, t, sol, status, path},
            ee = $ExampleBoundaryDefaults[key];
            entries = ee[[1]]; exits = ee[[2]];
            s = getExampleScenario[key, entries, exits];
            {t, sol} = AbsoluteTiming @ CheckAbort[
                Block[{$RecursionLimit = 16384, $IterationLimit = 16384},
                    TimeConstrained[Quiet @ solveScenario[s], $timeoutSec, $TimedOut]],
                $Aborted];
            status = Which[
                sol === $TimedOut, "TIMEOUT",
                sol === $Aborted, "ABORTED",
                Head[sol] === TerminatedEvaluation, "TERMINATED",
                True, "OK"
            ];
            path = If[status === "OK",
                storeSolution[key, entries, exits, sol],
                None
            ];
            AppendTo[results, <|"Key" -> key, "Status" -> status,
                "Time" -> t, "Path" -> path|>];
            Print[
                StringPadRight[ToString[key], 38], " ",
                StringPadRight[status, 12], " ",
                StringPadRight[ToString[NumberForm[N[t], {6, 2}]] <> "s", 12],
                If[path =!= None, " -> " <> FileBaseName[path], ""]
            ];
        ]
    ],
    $selectedKeys
];

okCount = Count[results, KeyValuePattern["Status" -> "OK"]];
skipCount = Length[results] - okCount;
Print[""];
Print["Wrote ", okCount, " cache file(s); skipped ", skipCount, " (timeout/abort)."];
