#!/usr/bin/env wolframscript
(* BenchmarkHelpers.wls -- Helper functions for the MFGraphs benchmark suite *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

(* --- Safe execution wrapper --- *)

SafeExecute::usage =
"SafeExecute[expr, timeout] evaluates expr with timeout protection and error handling.
Returns <|\"Time\" -> seconds, \"Memory\" -> bytes, \"Status\" -> string, \"Result\" -> value|>.
Status is \"OK\", \"TIMEOUT\", or \"FAILED\".";

SetAttributes[SafeExecute, HoldFirst];

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
solKey should be \"AssoCritical\", \"AssoNonCritical\", or \"AssoMonotone\".";

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

(* --- Solution hashing and verification --- *)

$ReferenceHashesPath = FileNameJoin[{Directory[], "MFGraphs", "Tests", "reference_solutions.json"}];

$ReferenceHashes = If[FileExistsQ[$ReferenceHashesPath],
  Quiet @ Check[
    Module[{data},
      data = Import[$ReferenceHashesPath, "RawJSON"];
      Association[
        (ToString[#["case"]] -> #) & /@ Lookup[data, "cases", {}]
      ]
    ],
    <||>
  ],
  <||>
];

HashSolution[asso_Association] :=
  Module[{pairs},
    pairs = SortBy[
      KeyValueMap[
        {ToString[#1, InputForm], If[NumericQ[#2], Round[N[#2], 10^-10], ToString[#2, InputForm]]} &,
        asso
      ],
      First
    ];
    Hash[pairs, "SHA256", "HexString"]
  ];

HashSolution[_] := None;

VerifySolution[caseKey_, asso_] :=
  Module[{ref, currentHash},
    ref = Lookup[$ReferenceHashes, ToString[caseKey], None];
    If[ref === None, Return["NO_REF"]];
    If[ref["solution_hash"] === None, Return["REF_NO_HASH"]];
    currentHash = HashSolution[asso];
    If[currentHash === None, Return["NO_SOLUTION"]];
    If[currentHash === ref["solution_hash"], "MATCH", "MISMATCH"]
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
  "HRF Scenario 1" -> {I1 -> 100, I2 -> 100, U1 -> 0, U2 -> 0},
  "Paper example" -> {I1 -> 80, I2 -> 40, S1 -> 2, S2 -> 1, S3 -> 3, S4 -> 1},
  13 -> {},
  "Braess split" -> {},
  "Braess congest" -> {},
  "New Braess" -> {},
  "Big Braess split" -> {},
  "Big Braess congest" -> {},
  "Grid0303" -> {},
  "Grid0404" -> {},
  "Grid0505" -> {},
  "Grid0707" -> {},
  "Grid0710" -> {},
  "Grid1010" -> {},
  "Grid1020" -> {}
|>;

GetParams[key_] := Lookup[$CaseParams, key, $DefaultParams];

(* --- Benchmark case tiers --- *)

$BenchmarkTiers = <|
  "small" -> {1, 2, 3, 4, 5, 6, 27},
  "core" -> {9, 10, 12, 14, 18},
  "stress" -> {7, 8, 104, "triangle with two exits"},
  "paper" -> {"HRF Scenario 1"},
  "inconsistent-switching" -> {11, 15, 17, "Paper example"},
  "large" -> {13, 19, 20, 21, 22, 23, "Jamaratv9",
              "Braess split", "Braess congest", "New Braess",
              "Big Braess split", "Big Braess congest",
              "Grid0303", "Grid0404", "Grid0505"},
  "vlarge" -> {"Grid0707", "Grid0710", "Grid1010", "Grid1020"}
|>;

$BenchmarkTierAliases = <|
  "medium" -> "core",
  "inconsistent" -> "inconsistent-switching"
|>;

$DefaultBenchmarkTiers = {"small", "core", "large", "vlarge"};

$TierTimeouts = <|
  "small" -> 60,
  "core" -> 300,
  "stress" -> 1800,
  "paper" -> 1800,
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

JSONSafeValue[value_] := Which[
  MatchQ[value, _Missing], Null,
  AssociationQ[value], Map[JSONSafeValue, value],
  ListQ[value], JSONSafeValue /@ value,
  NumberQ[value] || StringQ[value] || BooleanQ[value] || value === Null, value,
  True, ToString[Unevaluated[value], InputForm]
];

ExportResultsJSON[results_List, path_String] :=
  Module[{cleanResults},
    cleanResults = JSONSafeValue /@ results;
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
      {tier, {"small", "core", "stress", "paper", "inconsistent-switching", "large", "vlarge"}}
    ];
  ];
