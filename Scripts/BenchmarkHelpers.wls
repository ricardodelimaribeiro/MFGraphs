#!/usr/bin/env wolframscript
(* BenchmarkHelpers.wls -- Helper functions for the MFGraphs benchmark suite *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

(* --- Safe execution wrapper --- *)

SafeExecute::usage =
"SafeExecute[expr, timeout] evaluates expr with timeout protection and error handling.
Returns <|\"Time\" -> seconds, \"Memory\" -> bytes, \"Status\" -> string, \"Result\" -> value|>.
Status is \"OK\", \"TIMEOUT\", or \"FAILED\".";

SafeExecute[expr_, timeout_:300] :=
  Module[{result, time, memBefore, memAfter, status},
    memBefore = MemoryInUse[];
    {time, result} = AbsoluteTiming[
      TimeConstrained[
        Check[expr, $Failed],
        timeout,
        $TimedOut
      ]
    ];
    memAfter = MaxMemoryUsed[];
    status = Which[
      result === $TimedOut, "TIMEOUT",
      result === $Failed, "FAILED",
      result === Null, "FAILED",
      True, "OK"
    ];
    <|
      "Time" -> time,
      "Memory" -> memAfter - memBefore,
      "Status" -> status,
      "Result" -> result
    |>
  ];

(* --- Residual computation --- *)

ComputeResidual::usage =
"ComputeResidual[eqs, solKey] computes the infinity-norm residual from the solver output.
solKey should be \"AssoCritical\" or \"AssoNonCritical\".";

ComputeResidual[eqs_Association, solKey_String:"AssoNonCritical"] :=
  Module[{sol, nlhs, nrhs},
    sol = Lookup[eqs, solKey, <||>];
    nlhs = Lookup[eqs, "Nlhs", {}];
    nrhs = Lookup[eqs, "Nrhs", {}];
    If[sol === <||> || nlhs === {} || nrhs === {},
      Missing["NoData"],
      Quiet @ Check[
        Norm[N[nlhs /. sol] - (nrhs /. sol), Infinity],
        Missing["ComputeError"]
      ]
    ]
  ];

(* --- Graph metadata extraction --- *)

GraphMetadata::usage =
"GraphMetadata[data] returns an Association with vertex count, edge count, entry/exit counts, and switching flag.";

GraphMetadata[data_Association] :=
  Module[{vl, am, ev, xv, sc},
    vl = Lookup[data, "Vertices List", {}];
    am = Lookup[data, "Adjacency Matrix", {}];
    ev = Lookup[data, "Entrance Vertices and Flows", {}];
    xv = Lookup[data, "Exit Vertices and Terminal Costs", {}];
    sc = Lookup[data, "Switching Costs", {}];
    <|
      "NumVertices" -> Length[vl],
      "NumEdges" -> If[am === {}, 0, Total[Flatten[am]]],
      "NumEntries" -> Length[ev],
      "NumExits" -> Length[xv],
      "HasSwitching" -> (sc =!= {})
    |>
  ];

(* --- Default parameter substitution rules --- *)

$DefaultParams = {
  I1 -> 100, I2 -> 50, I3 -> 30,
  U1 -> 0, U2 -> 0, U3 -> 0,
  S1 -> 1, S2 -> 1, S3 -> 1, S4 -> 1,
  S5 -> 1, S6 -> 1, S7 -> 1, S8 -> 1,
  S9 -> 1, S10 -> 1, S11 -> 1, S12 -> 1,
  S13 -> 1, S14 -> 1, S15 -> 1, S16 -> 1
};

(* Per-case overrides: empty list means already numeric *)
$CaseParams = <|
  7 -> {I1 -> 80, U1 -> 0, U2 -> 10},
  8 -> {I1 -> 100, U1 -> 0, U2 -> 10, S1 -> 2, S2 -> 3, S3 -> 2, S4 -> 1, S5 -> 3, S6 -> 1},
  9 -> {I1 -> 100, I2 -> 40, U1 -> 0},
  10 -> {I1 -> 100, I2 -> 40, U1 -> 0, S1 -> 2, S2 -> 3, S3 -> 2, S4 -> 1, S5 -> 3, S6 -> 1},
  12 -> {I1 -> 80, U1 -> 0},
  14 -> {I1 -> 80, U1 -> 0, S1 -> 2, S2 -> 1},
  15 -> {I1 -> 80, I2 -> 40, U1 -> 0, S1 -> 2, S2 -> 1},
  17 -> {I1 -> 80, U1 -> 0, U2 -> 10, S1 -> 2, S2 -> 1},
  18 -> {I1 -> 80, U1 -> 0, S1 -> 2, S2 -> 1},
  104 -> {I1 -> 80, I2 -> 40},
  "triangle with two exits" -> {I1 -> 80, U1 -> 0, U2 -> 10},
  "Paper example" -> {I1 -> 80, I2 -> 40, S1 -> 2, S2 -> 1, S3 -> 3, S4 -> 1},
  13 -> {},
  "Braess split" -> {},
  "Braess congest" -> {},
  "New Braess" -> {},
  "Big Braess split" -> {},
  "Big Braess congest" -> {},
  "Grid0303" -> {},
  "Grid1020" -> {}
|>;

GetParams[key_] := Lookup[$CaseParams, key, $DefaultParams];

(* --- Benchmark case tiers --- *)

$BenchmarkTiers = <|
  "small" -> {1, 2, 3, 4, 5, 6, 27},
  "core" -> {9, 10, 12, 14, 15, 18, "Paper example"},
  "stress" -> {7, 8, 11, 17, 104, "triangle with two exits"},
  "large" -> {13, 19, 20, 21, 22, 23, "Jamaratv9",
              "Braess split", "Braess congest", "New Braess",
              "Big Braess split", "Big Braess congest", "Grid0303"},
  "vlarge" -> {"Grid1020"}
|>;

$BenchmarkTierAliases = <|
  "medium" -> "core"
|>;

$DefaultBenchmarkTiers = {"small", "core", "large", "vlarge"};

$TierTimeouts = <|
  "small" -> 60,
  "core" -> 300,
  "stress" -> 1800,
  "large" -> 900,
  "vlarge" -> 1800
|>;

(* --- Results formatting --- *)

FormatTime[t_?NumericQ] := ToString[CForm[N[t]]];
FormatTime[_] := "N/A";

FormatMemoryMB[bytes_?NumericQ] := ToString[CForm[N[bytes / 1024.^2]]];
FormatMemoryMB[_] := "N/A";

(* --- Export helpers --- *)

ExportResultsCSV[results_List, path_String] :=
  Module[{ds},
    ds = Dataset[results];
    Export[path, ds, "CSV"];
    Print["Results exported to ", path];
  ];

ExportResultsJSON[results_List, path_String] :=
  Module[{cleanResults},
    cleanResults = results /. _Missing -> Null;
    Export[path, cleanResults, "JSON"];
    Print["Results exported to ", path];
  ];

PrintResultsSummary[results_List] :=
  Module[{grouped, stats},
    Print["\n", Style["=== Benchmark Results Summary ===", Bold]];
    Print["Total cases run: ", Length[results]];
    grouped = GroupBy[results, #["Status"] &];
    Print["  OK:      ", Length @ Lookup[grouped, "OK", {}]];
    Print["  FAILED:  ", Length @ Lookup[grouped, "FAILED", {}]];
    Print["  TIMEOUT: ", Length @ Lookup[grouped, "TIMEOUT", {}]];
    Print["  SKIPPED: ", Length @ Lookup[grouped, "SKIPPED", {}]];
    Print[];
    (* Print per-tier summary *)
    Do[
      stats = Select[results, #["Tier"] === tier &];
      If[Length[stats] > 0,
        Print[Style["--- Tier: " <> tier <> " ---", Bold]];
        Do[
          Print["  ", r["Key"], " | ", r["Solver"], " | ",
                FormatTime[r["SolveTime"]], "s | ",
                FormatMemoryMB[r["SolveMemory"]], " MB | ",
                r["Status"]],
          {r, stats}
        ];
        Print[];
      ],
      {tier, {"small", "core", "stress", "large", "vlarge"}}
    ];
  ];
