#!/usr/bin/env wolframscript

PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];

jsonPath = FileNameJoin[{Directory[], "Results", "benchmark_latest.json"}];
tiersToRun = {"small", "core"};
results = Flatten @ Reap[
    Scan[
        Function[tier,
            Module[{run, tierRows},
                Print["[CI] Running tier: ", tier];
                run = RunProcess[{
                    "wolframscript",
                    "-file", "Scripts/BenchmarkSuite.wls",
                    tier,
                    "isolation=inprocess",
                    "timeout=10"
                }];
                If[Lookup[run, "ExitCode", 1] =!= 0,
                    Print["[FAIL] BenchmarkSuite tier ", tier, " failed with exit code ", Lookup[run, "ExitCode", Missing["Unknown"]], "."];
                    If[StringLength @ Lookup[run, "StandardError", ""] > 0,
                        Print["[stderr] ", Lookup[run, "StandardError", ""]]
                    ];
                    Exit[1];
                ];

                If[!FileExistsQ[jsonPath],
                    Print["[FAIL] Benchmark JSON not found at ", jsonPath, " after tier ", tier, "."];
                    Exit[1];
                ];

                tierRows = Quiet @ Check[Import[jsonPath, "RawJSON"], $Failed];
                If[tierRows === $Failed || !ListQ[tierRows],
                    Print["[FAIL] Could not parse benchmark JSON after tier ", tier, "."];
                    Exit[1];
                ];
                Sow /@ tierRows;
            ]
        ],
        tiersToRun
    ]
][[2, 1]];

failures = {};
AppendFailure[msg_] := AppendTo[failures, msg];

Scan[
    Function[row,
        Module[{key, solver, status, resultKind, feasibility, message,
                hasResultKind, hasFeasibility, hasMessage},
            If[!AssociationQ[row],
                AppendFailure["Non-association benchmark row encountered."];
                Return[];
            ];

            key = ToString @ Lookup[row, "Key", Missing["UnknownKey"]];
            solver = ToString @ Lookup[row, "Solver", Missing["UnknownSolver"]];
            status = Lookup[row, "Status", Missing["UnknownStatus"]];

            hasResultKind = KeyExistsQ[row, "ResultKind"];
            hasFeasibility = KeyExistsQ[row, "Feasibility"];
            hasMessage = KeyExistsQ[row, "Message"];

            resultKind = Lookup[row, "ResultKind", Missing["MissingResultKind"]];
            feasibility = Lookup[row, "Feasibility", Missing["MissingFeasibility"]];
            message = Lookup[row, "Message", Missing["MissingMessage"]];

            (* Rule 0: Benchmark schema must expose truth fields at top-level. *)
            If[!(hasResultKind && hasFeasibility && hasMessage),
                AppendFailure[
                    "[FAIL] Missing truth-schema fields on row " <> key <> " (" <> solver <> ")."
                ];
            ];

            (* Rule 1: No kernel crashes. *)
            If[status === "FAILED",
                AppendFailure[
                    "[FAIL] Kernel crash detected on row " <> key <> " (" <> solver <> ")."
                ];
            ];

            (* Rule 2: Successful execution must include materialized truth status. *)
            If[status === "OK" && (resultKind === Null || MatchQ[resultKind, _Missing]),
                AppendFailure[
                    "[FAIL] Missing ResultKind for successful row " <> key <> " (" <> solver <> ")."
                ];
            ];
            If[status === "OK" && (feasibility === Null || MatchQ[feasibility, _Missing]),
                AppendFailure[
                    "[FAIL] Missing Feasibility for successful row " <> key <> " (" <> solver <> ")."
                ];
            ];
            If[status === "OK" && MatchQ[message, _Missing],
                AppendFailure[
                    "[FAIL] Missing Message for successful row " <> key <> " (" <> solver <> ")."
                ];
            ];

            (* Rule 3: Monotone must not fake success on strict residual violation. *)
            If[
                solver === "Monotone" && status === "OK" &&
                message === "ResidualExceedsTolerance" &&
                resultKind === "Success",
                AppendFailure[
                    "[FAIL] Truth-gate violation on row " <> key <> " (" <> solver <> "): ResidualExceedsTolerance labeled Success."
                ];
            ];
        ]
    ],
    results
];

If[Length[failures] > 0,
    Print["[CI] FAILURES DETECTED: ", Length[failures]];
    Scan[Print, failures];
    Exit[1],
    Print["[CI] PASS: benchmark truth schema checks succeeded for ", Length[results], " rows."];
    Exit[0]
];
