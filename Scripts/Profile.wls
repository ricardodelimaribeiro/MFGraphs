#!/usr/bin/env wolframscript
(* Profile.wls -- Run a single benchmark case with profiling enabled *)

Print["MFGraphs Profiler"];
Print["================="];
Print["Date: ", DateString[]];
Print[];

(* --- Setup package path --- *)
$RepoRoot = "/Users/ribeirrd/Documents/GitHub/MFGraphs";
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

(* --- Load helpers and profiler --- *)
Get[FileNameJoin[{$RepoRoot, "Scripts", "BenchmarkHelpers.wls"}]];
Get[FileNameJoin[{$RepoRoot, "Scripts", "ProfileInstrument.wls"}]];

(* --- Configuration --- *)
$CaseToProfile = 7;
$SolverToProfile = "CriticalCongestion";

(* --- Main execution --- *)
Module[{tier, timeout, data, d2eResult, d2e, solverResult, report},
  tier = GetTier[$CaseToProfile];
  timeout = 900;

  Print["--- Profiling Case: ", $CaseToProfile, " (tier: ", tier, ") ---"];
  Print["Solver: ", $SolverToProfile];

  (* Load data *)
  data = Quiet @ Check[GetExampleData[$CaseToProfile], $Failed];
  If[data === $Failed,
    Print["  FAILED to load GetExampleData[", $CaseToProfile, "]"];
    Quit[];
  ];

  (* Convert to equations *)
  Print["  Running DataToEquations..."];
  d2eResult = SafeExecute[DataToEquations[data], timeout];
  Print["    D2E Status: ", d2eResult["Status"]];

  If[d2eResult["Status"] =!= "OK",
    Print["  DataToEquations failed, cannot proceed with profiling."];
    Quit[];
  ];
  d2e = d2eResult["Result"];

  (* Run solver with profiling *)
  Print["  Running ", $SolverToProfile, " with profiling..."];
  
  InitProfile[];
  InstallProfileWrappers[];

  Print["Entering profiled CriticalCongestionSolver block..."];
  solverResult = SafeExecute[CriticalCongestionSolver[d2e], timeout];
  Print["...left profiled CriticalCongestionSolver block."];
  
  UninstallProfileWrappers[];

  Print["    Solver Status: ", solverResult["Status"]];
  Print["    Total Time: ", FormatTime[solverResult["Time"]], "s"];
  Print[];

  (* Generate and print report *)
  Print[Style["--- Profiling Report ---", Bold]];
  report = ProfileReportMarkdown[];
  Print[report];
];

Print["
Profiling complete at ", DateString[]];
