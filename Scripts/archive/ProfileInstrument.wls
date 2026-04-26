#!/usr/bin/env wolframscript
(* ProfileInstrument.wls -- Non-invasive profiling instrumentation for MFGraphs *)
(* Load AFTER Needs["MFGraphs`"] *)

(* --- Profile accumulators --- *)

$ProfileTimings = <||>;
$ProfileCounts = <||>;
$ProfileExtras = <||>;

InitProfile[] := (
  $ProfileTimings = <|
    "DNFReduce" -> 0., "DNFReduce_Solve" -> 0., "DNFReduce_Reduce" -> 0.,
    "DNFReduce_Simplify" -> 0., "TripleStep" -> 0., "TripleClean" -> 0.,
    "DNFSolveStep" -> 0.,
    "GradientProjection" -> 0., "PseudoInverse" -> 0.,
    "MFGSystemSolver" -> 0., "MFGPreprocessing" -> 0.
  |>;
  $ProfileCounts = <|
    "DNFReduce" -> 0, "DNFReduce_Solve" -> 0, "DNFReduce_Reduce" -> 0,
    "DNFReduce_Simplify" -> 0, "TripleStep" -> 0, "TripleClean" -> 0,
    "DNFSolveStep" -> 0,
    "GradientProjection" -> 0, "PseudoInverse" -> 0,
    "MFGSystemSolver" -> 0, "MFGPreprocessing" -> 0
  |>;
  $ProfileExtras = <|
    "DNFReduce_MaxDepth" -> 0,
    "DNFReduce_CurrentDepth" -> 0,
    "TripleClean_MaxIters" -> 0
  |>;
);

(* --- Profile report generation --- *)

ProfileReport[] :=
  Module[{rows},
    rows = KeyValueMap[
      Function[{key, time},
        <|
          "Function" -> key,
          "TotalTime" -> time,
          "CallCount" -> Lookup[$ProfileCounts, key, 0],
          "AvgTime" -> If[$ProfileCounts[key] > 0, time / $ProfileCounts[key], 0.]
        |>
      ],
      $ProfileTimings
    ];
    (* Add extras *)
    rows = Append[rows, <|
      "Function" -> "DNFReduce_MaxDepth",
      "TotalTime" -> Missing["N/A"],
      "CallCount" -> $ProfileExtras["DNFReduce_MaxDepth"],
      "AvgTime" -> Missing["N/A"]
    |>];
    rows = Append[rows, <|
      "Function" -> "TripleClean_MaxIters",
      "TotalTime" -> Missing["N/A"],
      "CallCount" -> $ProfileExtras["TripleClean_MaxIters"],
      "AvgTime" -> Missing["N/A"]
    |>];
    Dataset[rows]
  ];

ProfileReportMarkdown[] :=
  Module[{rows, header, lines},
    rows = Normal @ ProfileReport[];
    header = "| Function | Total Time (s) | Call Count | Avg Time (s) |";
    lines = {header, "|---|---|---|---|"};
    Do[
      AppendTo[lines,
        "| " <> ToString[r["Function"]] <> " | " <>
        ToString[NumberForm[r["TotalTime"], {6, 4}]] <> " | " <>
        ToString[r["CallCount"]] <> " | " <>
        ToString[NumberForm[r["AvgTime"], {6, 6}]] <> " |"
      ],
      {r, rows}
    ];
    StringRiffle[lines, "\n"]
  ];

(* --- Instrumentation via DownValues wrapping --- *)

(* Flags to prevent recursive instrumentation *)
$InProfileDNFReduce = False;
$InProfileTripleStep = False;
$InProfileTripleClean = False;
$InProfileDNFSolveStep = False;
$InProfileGradProj = False;
$InProfileMFGSystemSolver = False;

InstallProfileWrappers[] :=
  Module[{},
    Print["Installing profiling wrappers..."];

    (* --- Wrap DNFReduce --- *)
    $OriginalDNFReduceDV = DownValues[DNFReduce];
    DownValues[DNFReduce] = Join[
      {
        HoldPattern[DNFReduce[xp_, arg_]] :>
          Module[{res, t},
            $ProfileCounts["DNFReduce"]++;
            $ProfileExtras["DNFReduce_CurrentDepth"]++;
            $ProfileExtras["DNFReduce_MaxDepth"] = Max[
              $ProfileExtras["DNFReduce_MaxDepth"],
              $ProfileExtras["DNFReduce_CurrentDepth"]
            ];
            t = AbsoluteTime[];
            Block[{$InProfileDNFReduce = True},
              With[{savedDV = DownValues[DNFReduce]},
                DownValues[DNFReduce] = $OriginalDNFReduceDV;
                res = DNFReduce[xp, arg];
                DownValues[DNFReduce] = savedDV;
              ]
            ];
            $ProfileTimings["DNFReduce"] += AbsoluteTime[] - t;
            $ProfileExtras["DNFReduce_CurrentDepth"]--;
            res
          ] /; !TrueQ[$InProfileDNFReduce]
      },
      $OriginalDNFReduceDV
    ];

    (* --- Wrap TripleStep --- *)
    $OriginalTripleStepDV = DownValues[TripleStep];
    DownValues[TripleStep] = Join[
      {
        HoldPattern[TripleStep[arg_]] :>
          Module[{res, t},
            $ProfileCounts["TripleStep"]++;
            t = AbsoluteTime[];
            Block[{$InProfileTripleStep = True},
              With[{savedDV = DownValues[TripleStep]},
                DownValues[TripleStep] = $OriginalTripleStepDV;
                res = TripleStep[arg];
                DownValues[TripleStep] = savedDV;
              ]
            ];
            $ProfileTimings["TripleStep"] += AbsoluteTime[] - t;
            res
          ] /; !TrueQ[$InProfileTripleStep]
      },
      $OriginalTripleStepDV
    ];

    (* --- Wrap TripleClean --- *)
    $OriginalTripleCleanDV = DownValues[TripleClean];
    DownValues[TripleClean] = Join[
      {
        HoldPattern[TripleClean[arg_]] :>
          Module[{res, t},
            $ProfileCounts["TripleClean"]++;
            t = AbsoluteTime[];
            Block[{$InProfileTripleClean = True},
              With[{savedDV = DownValues[TripleClean]},
                DownValues[TripleClean] = $OriginalTripleCleanDV;
                res = TripleClean[arg];
                DownValues[TripleClean] = savedDV;
              ]
            ];
            $ProfileTimings["TripleClean"] += AbsoluteTime[] - t;
            res
          ] /; !TrueQ[$InProfileTripleClean]
      },
      $OriginalTripleCleanDV
    ];

    (* --- Wrap DNFSolveStep --- *)
    $OriginalDNFSolveStepDV = DownValues[DNFSolveStep];
    DownValues[DNFSolveStep] = Join[
      {
        HoldPattern[DNFSolveStep[arg_]] :>
          Module[{res, t},
            $ProfileCounts["DNFSolveStep"]++;
            t = AbsoluteTime[];
            Block[{$InProfileDNFSolveStep = True},
              With[{savedDV = DownValues[DNFSolveStep]},
                DownValues[DNFSolveStep] = $OriginalDNFSolveStepDV;
                res = DNFSolveStep[arg];
                DownValues[DNFSolveStep] = savedDV;
              ]
            ];
            $ProfileTimings["DNFSolveStep"] += AbsoluteTime[] - t;
            res
          ] /; !TrueQ[$InProfileDNFSolveStep]
      },
      $OriginalDNFSolveStepDV
    ];

    (* --- Wrap GradientProjection --- *)
    $OriginalGradProjDV = DownValues[GradientProjection];
    DownValues[GradientProjection] = Join[
      {
        HoldPattern[GradientProjection[x_?NumberVectorQ, A_, dim_, At_]] :>
          Module[{res, t},
            $ProfileCounts["GradientProjection"]++;
            t = AbsoluteTime[];
            Block[{$InProfileGradProj = True},
              With[{savedDV = DownValues[GradientProjection]},
                DownValues[GradientProjection] = $OriginalGradProjDV;
                res = GradientProjection[x, A, dim, At];
                DownValues[GradientProjection] = savedDV;
              ]
            ];
            $ProfileTimings["GradientProjection"] += AbsoluteTime[] - t;
            res
          ] /; !TrueQ[$InProfileGradProj]
      },
      $OriginalGradProjDV
    ];

    (* --- Wrap MFGSystemSolver --- *)
    $OriginalMFGSystemSolverDV = DownValues[MFGSystemSolver];
    DownValues[MFGSystemSolver] = Join[
        {
            HoldPattern[MFGSystemSolver[Eqs_][approxJs_]] :>
                Module[{res, t},
                    $ProfileCounts["MFGSystemSolver"]++;
                    t = AbsoluteTime[];
                    Block[{$InProfileMFGSystemSolver = True},
                        With[{savedDV = DownValues[MFGSystemSolver]},
                            DownValues[MFGSystemSolver] = $OriginalMFGSystemSolverDV;
                            res = MFGSystemSolver[Eqs][approxJs];
                            DownValues[MFGSystemSolver] = savedDV;
                        ]
                    ];
                    $ProfileTimings["MFGSystemSolver"] += AbsoluteTime[] - t;
                    res
                ] /; !TrueQ[$InProfileMFGSystemSolver]
        },
        $OriginalMFGSystemSolverDV
    ];

    Print["Profiling wrappers installed."];
  ];

UninstallProfileWrappers[] :=
  Module[{},
    If[ValueQ[$OriginalDNFReduceDV], DownValues[DNFReduce] = $OriginalDNFReduceDV];
    If[ValueQ[$OriginalTripleStepDV], DownValues[TripleStep] = $OriginalTripleStepDV];
    If[ValueQ[$OriginalTripleCleanDV], DownValues[TripleClean] = $OriginalTripleCleanDV];
    If[ValueQ[$OriginalDNFSolveStepDV], DownValues[DNFSolveStep] = $OriginalDNFSolveStepDV];
    If[ValueQ[$OriginalGradProjDV], DownValues[GradientProjection] = $OriginalGradProjDV];
    If[ValueQ[$OriginalMFGSystemSolverDV], DownValues[MFGSystemSolver] = $OriginalMFGSystemSolverDV];
    Print["Profiling wrappers removed."];
  ];

(* --- Profiled solver convenience wrappers --- *)

ProfiledCriticalCongestion[d2e_] :=
  Module[{result, t0},
    InitProfile[];
    InstallProfileWrappers[];
    t0 = AbsoluteTime[];
    result = CriticalCongestionSolver[d2e];
    $ProfileTimings["CriticalCongestionSolver_Total"] = AbsoluteTime[] - t0;
    UninstallProfileWrappers[];
    result
  ];

Print["ProfileInstrument.wls loaded successfully."];
Print["Usage: InitProfile[]; InstallProfileWrappers[]; <run solver>; ProfileReport[]"];
