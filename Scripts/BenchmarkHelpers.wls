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

SafeExecuteIsolated::usage =
"SafeExecuteIsolated[expr, timeout] evaluates expr in an isolated child wolframscript process.
Returns the same status envelope as SafeExecute and prevents child-kernel crashes from
terminating the benchmark master process.";

ResolveWolframScriptKernelPath[] :=
  Module[{parentExecutable},
    parentExecutable = Quiet @ Check[First[$CommandLine], $Failed];
    Which[
      StringQ[parentExecutable] &&
        StringContainsQ[ToLowerCase[parentExecutable], "wolframscript"] &&
        FileExistsQ[parentExecutable],
        parentExecutable,
      FileExistsQ["/usr/local/bin/wolframscript"],
        "/usr/local/bin/wolframscript",
      True,
        "wolframscript"
    ]
  ];

Options[SafeExecuteIsolated] = {
  "KernelPath" -> ResolveWolframScriptKernelPath[],
  "ProcessTimeBuffer" -> 5
};

$SafeExecuteMode = "InProcess";
$SafeExecuteKernelPath = Automatic;

SetAttributes[SafeExecuteIsolated, HoldFirst];
SetAttributes[SafeExecuteInProcess, HoldFirst];

SafeExecuteInProcess[expr_, timeout_:300] :=
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

SafeExecuteIsolated[expr_, timeout_:300, opts:OptionsPattern[]] :=
  Module[{kernelPath, buffer, scriptFile, outFile, runResult = <||>, runTime = 0.,
      rawResult = Missing["NotAvailable"], status = "FAILED", result = Missing["Failed"],
      stderr = "", heldExpr, concretizedExpr, exprString, scriptLines, processTimeout,
      proc, startTime, timedOut = False, info = <||>, stdout = "", startOK = False,
      ignored, procStatus = "Unknown", rawPayload = Missing["NotAvailable"],
      childTime = Missing["NotAvailable"], childMem = Missing["NotAvailable"]},

    kernelPath = OptionValue["KernelPath"];
    buffer = OptionValue["ProcessTimeBuffer"];
    processTimeout = Max[1, timeout + buffer];

    scriptFile = CreateFile[];
    outFile = CreateFile[];

    heldExpr = HoldComplete[expr];
    (* Inline caller-bound symbols into the child script. Two pitfalls this avoids:
       1. Evaluate[s] inside HoldComplete is inert — must use RuleCondition to force RHS eval.
       2. Symbols can land outside Global`: Module creates local contexts (e.g. foo$123`),
          and name collisions with MFGraphs` exports put user variables in MFGraphs`.
       We inline any non-System symbol that has OwnValues. *)
    concretizedExpr = heldExpr /. s_Symbol /;
        (!StringStartsQ[Context[s], "System`"] && OwnValues[s] =!= {}) :>
      RuleCondition[s /. OwnValues[s], True];
    exprString = StringReplace[
      ToString[concretizedExpr, InputForm],
      {
        StartOfString ~~ "HoldComplete[" -> "",
        "]" ~~ EndOfString -> ""
      }
    ];

    scriptLines = {
      "SetDirectory[" <> ToString[InputForm[Directory[]]] <> "];",
      "PrependTo[$Path, Directory[]];",
      "Needs[\"MFGraphs`\"];",
      "Quiet[Check[ResetMaxMemoryUsed[], Null]];",
      "{innerTime, rawRes} = AbsoluteTiming[",
      "  TimeConstrained[",
      "    Check[",
      "      " <> exprString <> ",",
      "      $Failed",
      "    ],",
      "    " <> ToString[InputForm[N[timeout]]] <> ",",
      "    \"TIMEOUT\"",
      "  ]",
      "];",
      "childMem = Quiet @ Check[MaxMemoryUsed[], Missing[\"NotAvailable\"]];",
      "payload = <|\"Result\" -> rawRes, \"InnerSolveTime\" -> innerTime, \"ChildMaxMemoryUsed\" -> childMem|>;",
      "Export[" <> ToString[InputForm[outFile]] <> ", payload, \"WXF\"];"
    };
    Export[scriptFile, StringRiffle[scriptLines, "\n"], "Text"];

    {runTime, ignored} = AbsoluteTiming[
      Quiet @ Check[
        proc = StartProcess[{kernelPath, "-file", scriptFile}];
        startOK = True,
        startOK = False
      ];
      If[startOK,
        startTime = AbsoluteTime[];
        procStatus = Quiet @ Check[ProcessStatus[proc], "Unknown"];
        While[
          procStatus === "Running" && AbsoluteTime[] - startTime < processTimeout,
          Pause[0.05]
          ;
          procStatus = Quiet @ Check[ProcessStatus[proc], "Unknown"]
        ];
        If[procStatus === "Running",
          timedOut = True;
          Quiet @ Check[KillProcess[proc], Null];
          Pause[0.1];
          procStatus = Quiet @ Check[ProcessStatus[proc], "Unknown"]
        ];
        info = Quiet @ Check[ProcessInformation[proc], <||>];
        If[procStatus =!= "Running",
          stdout = Quiet @ Check[ReadString[proc["StandardOutput"]], ""];
          stderr = Quiet @ Check[ReadString[proc["StandardError"]], ""];
        ];
      ];
      runResult = <|
        "ExitCode" -> Lookup[info, "ExitCode", If[startOK, Missing["NotAvailable"], 1]],
        "StandardOutput" -> stdout,
        "StandardError" -> stderr
      |>;
    ];

    If[timedOut ||
       Lookup[runResult, "ExitCode", 1] =!= 0 ||
       !FileExistsQ[outFile] ||
       Quiet @ Check[FileByteCount[outFile], 0] == 0,
      status =
        If[timedOut, "TIMEOUT", "FAILED"];
      result = Missing[If[status === "TIMEOUT", "Timeout", "KernelPanic"]],
      rawPayload = Quiet @ Check[Import[outFile, "WXF"], $Failed];
      If[AssociationQ[rawPayload] && KeyExistsQ[rawPayload, "Result"],
        rawResult = Lookup[rawPayload, "Result", $Failed];
        childTime = Lookup[rawPayload, "InnerSolveTime", Missing["NotAvailable"]];
        childMem = Lookup[rawPayload, "ChildMaxMemoryUsed", Missing["NotAvailable"]];
        ,
        rawResult = $Failed;
        childTime = Missing["NotAvailable"];
        childMem = Missing["NotAvailable"]
      ];
      Which[
        rawResult === "TIMEOUT",
          status = "TIMEOUT";
          result = Missing["Timeout"],
        rawResult === $Failed || rawResult === $Aborted || rawResult === Null,
          status = "FAILED";
          result = Missing["Failed"],
        True,
          status = "OK";
          result = rawResult
      ]
    ];

    Quiet @ Check[DeleteFile[scriptFile], Null];
    Quiet @ Check[DeleteFile[outFile], Null];

    <|
      "Time" -> runTime,
      "WallTime" -> runTime,
      "Memory" -> childMem,
      "Status" -> status,
      "Result" -> result,
      "InnerSolveTime" -> childTime,
      "ChildMaxMemoryUsed" -> childMem,
      "ExecutionMode" -> "Process",
      "ExitCode" -> Lookup[runResult, "ExitCode", Missing["NotAvailable"]],
      "StandardError" -> stderr
    |>
  ];

SafeExecute[expr_, timeout_:300] :=
  Module[{kernelPathOption},
    kernelPathOption = Replace[$SafeExecuteKernelPath, {
      s_String /; StringLength[StringTrim[s]] > 0 :>
        Sequence["KernelPath" -> StringTrim[s]],
      _ :> Sequence[]
    }];
    If[$SafeExecuteMode === "Process",
      SafeExecuteIsolated[expr, timeout, kernelPathOption],
      SafeExecuteInProcess[expr, timeout]
    ]
  ];

(* --- Residual computation --- *)

ComputeResidual::usage =
"ComputeResidual[eqs, solKey] computes the infinity-norm residual from the solver output.
solKey should be \"AssoCritical\", \"AssoNonCritical\", or \"AssoMonotone\".";

ComputeResidual[eqs_Association, solKey_String:"AssoNonCritical"] :=
  Module[{sol, nlhs, nrhs, diff},
    sol = Lookup[eqs, solKey, <||>];
    nlhs = Lookup[eqs, "Nlhs", {}];
    nrhs = Lookup[eqs, "Nrhs", {}];
    If[
      !AssociationQ[sol] || sol === <||> ||
      !ListQ[nlhs] || !ListQ[nrhs] || nlhs === {} || nrhs === {} ||
      Length[nlhs] =!= Length[nrhs],
      Return[Missing["NoData"], Module]
    ];
    diff = Quiet @ Check[(nlhs - nrhs) /. sol, $Failed];
    If[diff === $Failed,
      Return[Missing["ComputeError"], Module]
    ];
    diff = Quiet @ Check[N @ diff, $Failed];
    If[diff === $Failed || !VectorQ[diff, NumericQ],
      Missing["NotComputable"],
      Quiet @ Check[N @ Norm[diff, Infinity], Missing["ComputeError"]]
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

(* --- DNF Benchmark Helpers --- *)

ExtractDNFInput::usage =
"ExtractDNFInput[eqs] extracts the boolean system (xp, sorted, initRules) that would be passed to DNFReduce.";

ExtractDNFInput[eqs_Association] :=
  Module[{preEqs, newSystem, initRules, sorted, xp},
    (* Run preprocessing if not already done *)
    preEqs = If[KeyExistsQ[eqs, "InitRules"],
      eqs,
      MFGPreprocessing[eqs]
    ];
    initRules = Lookup[preEqs, "InitRules", <||>];
    newSystem = Lookup[preEqs, "NewSystem", {True, True, True}];

    (* Simulate what MFGSystemSolver does before DNFSolveStep: *)
    (* Apply zero-flow (critical congestion) substitution *)
    Module[{js, approxJs, costpluscurrents, Ncpc},
      js = Lookup[preEqs, "js", {}];
      approxJs = AssociationThread[js, 0 js];
      costpluscurrents = Lookup[preEqs, "costpluscurrents", {}];
      Ncpc = Expand /@ (costpluscurrents /. approxJs);
      initRules = Expand /@ (initRules /. Ncpc);
      newSystem = newSystem /. Ncpc;
    ];

    (* TripleClean *)
    {newSystem, initRules} = TripleClean[{newSystem, initRules}];

    (* Build the arguments that DNFSolveStep passes to DNFReduce *)
    If[newSystem[[3]] === True,
      (* No disjunctions — nothing interesting to compare *)
      {True, True, initRules, "trivial"},
      sorted = DeduplicateByComplexity[newSystem[[3]]];
      xp = And @@ Most[newSystem];
      {xp, sorted, initRules, "nontrivial"}
    ]
  ];

(* Count the number of disjuncts in a DNF expression *)
CountDisjuncts[expr_Or] := Length[expr];
CountDisjuncts[False] := 0;
CountDisjuncts[True] := 1;
CountDisjuncts[_] := 1;

(* Count conjuncts per disjunct *)
ConjunctsPerDisjunct[expr_Or] := Length /@ (List @@ expr) /. {n_Integer :> n, x_ :> 1};
ConjunctsPerDisjunct[expr_And] := {Length[expr]};
ConjunctsPerDisjunct[_] := {1};

(* --- Solver Result Verification Utilities --- *)

SuccessResultQ::usage =
"SuccessResultQ[result] returns True if the solver result Association indicates a successful solve with a valid solution.";

SuccessResultQ[result_Association] :=
  Lookup[result, "ResultKind", "Failure"] === "Success" &&
  AssociationQ[Lookup[result, "Solution", Missing["NotAvailable"]]];

VectorDifference::usage =
"VectorDifference[a, b] returns the infinity-norm difference between two numeric vectors, or Infinity if not comparable.";

VectorDifference[a_, b_] /; VectorQ[a, NumericQ] && VectorQ[b, NumericQ] :=
  Max[Abs[a - b]];
VectorDifference[_, _] := Infinity;

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
  "chain with two exits" -> {I1 -> 80, U1 -> 0, U2 -> 10},
  105 -> {I1 -> 80, U1 -> 0, U2 -> 10},
  "HRF Scenario 1" -> {I1 -> 100, I2 -> 100, U1 -> 0, U2 -> 0},
  "Jamaratv9" -> {I1 -> 100, I2 -> 50, U1 -> 0, U2 -> 0, U3 -> 0},
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

(* Named cases for specific benchmark scripts (CompareDNF, etc.) *)
$NamedBenchmarkCases = <|
  "Y_noswitch" -> {7, {I1 -> 50, U1 -> 0, U2 -> 0}},
  "Y_switch" -> {8, {I1 -> 100, U1 -> 0, U2 -> 10, S1 -> 2, S2 -> 3, S3 -> 2, S4 -> 1, S5 -> 3, S6 -> 1}},
  "Attraction_noswitch" -> {12, {I1 -> 100, U1 -> 0}},
  "Attraction_switch" -> {11, {I1 -> 100, U1 -> 0, S1 -> 1, S2 -> 1, S3 -> 1, S4 -> 1}},
  "Triangle_switch" -> {14, {I1 -> 100, U1 -> 0, U2 -> 0, S1 -> 1, S2 -> 1, S3 -> 1, S4 -> 1, S5 -> 1, S6 -> 1}},
  "Braess" -> {"Braess congest", {}},
  "Paper" -> {"Paper example", {I1 -> 100, I2 -> 50, U1 -> 0, U2 -> 0, S1 -> 1, S2 -> 1, S3 -> 1, S4 -> 1, S5 -> 1, S6 -> 1, S7 -> 1, S8 -> 1, S9 -> 1, S10 -> 1, S11 -> 1, S12 -> 1, S13 -> 1, S14 -> 1, S15 -> 1, S16 -> 1}},
  "Multi_20" -> {20, {I1 -> 80, I2 -> 60, U1 -> 0, U2 -> 0, U3 -> 0}},
  "Grid0303" -> {"Grid0303", {}},
  "Grid0505" -> {"Grid0505", {}},
  "Grid0710" -> {"Grid0710", {}},
  "Grid1010" -> {"Grid1010", {}},
  "Grid1020" -> {"Grid1020", {}}
|>;

GetParams[key_] := Lookup[$CaseParams, key, $DefaultParams];

(* --- Benchmark case tiers --- *)

$BenchmarkTiers = <|
  "small" -> {1, 2, 3, 4, 5, 6, 27},
  "core" -> {9, 10, 12, 14, 18, "chain with two exits"},
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

(* --- Benchmark Profiles --- *)

$BenchmarkProfiles = <|
  "critical-core" -> <|
    "Tiers" -> {"small", "core"},
    "Solvers" -> {"critical"}
  |>,
  "nonlinear-kkt" -> <|
    "Tiers" -> {"core", "stress", "paper"},
    "Solvers" -> {"nonlinear"}
  |>,
  "dnf-paper" -> <|
    "Tiers" -> {"paper", "inconsistent-switching"},
    "Solvers" -> {"all"}
  |>
|>;

(* --- Results formatting --- *)

DisplaySolverName[r_Association] :=
  Module[{base, details = {}},
    base = Switch[r["Solver"],
      "CriticalCongestion", "Critical",
      "NonLinearSolver", "NonLinear",
      "Monotone", "Monotone",
      "ALL", "Combined",
      "DataToEquations", "Preprocessing",
      _, r["Solver"]
    ];

    (* Add seeding info *)
    If[TrueQ[r["UsedCriticalSeed"]], AppendTo[details, "seeded"]];

    (* Add critical strategy details if available *)
    If[r["Solver"] === "CriticalCongestion" && KeyExistsQ[r, "NumericBackendUsed"],
      AppendTo[details, If[TrueQ[r["NumericBackendUsed"]], "numeric", "symbolic"]]
    ];

    If[details === {},
      base,
      base <> " (" <> StringRiffle[details, ", "] <> ")"
    ]
  ];

DisplaySolverName[name_String] := DisplaySolverName[<|"Solver" -> name|>]

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
          Print["  ", r["Key"], " | ", DisplaySolverName[r], " | ",
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
