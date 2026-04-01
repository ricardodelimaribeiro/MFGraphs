#!/usr/bin/env wolframscript
(* GenerateReferenceHashes.wls — Produce MFGraphs/Tests/reference_solutions.json

   For each benchmark case, solve CriticalCongestionSolver and hash the solution
   so future runs can verify correctness without comparing full associations.

   Usage:
     wolframscript -file Scripts/GenerateReferenceHashes.wls
     wolframscript -file Scripts/GenerateReferenceHashes.wls small
     wolframscript -file Scripts/GenerateReferenceHashes.wls medium
*)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

SetDirectory[DirectoryName[$InputFileName] // ParentDirectory];

Needs["MFGraphs`"];
Get[FileNameJoin[{Directory[], "Scripts", "BenchmarkHelpers.wls"}]];
$MFGraphsVerbose = False;

(* --- Deterministic solution hashing --- *)

HashSolution[asso_Association] :=
  Module[{pairs},
    (* Sort by string key for determinism, round numerics to avoid float drift *)
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

HashParams[params_List] := Hash[Sort[params], "SHA256", "HexString"];
HashParams[_] := Hash[{}, "SHA256", "HexString"];

(* --- Select cases --- *)

tier = SelectFirst[Rest[$ScriptCommandLine], !StringStartsQ[#, "--"]&, "all"];
cases = If[tier === "all",
  Flatten[Values[$BenchmarkTiers]],
  Lookup[$BenchmarkTiers, tier, {}]
];

Print["Generating reference hashes for ", Length[cases], " cases (tier: ", tier, ")...\n"];

(* --- Solve and hash --- *)

entries = {};

Do[
  Module[{params, data, eqs, result, asso, solHash, paramsHash, status, nKeys},
    params = GetParams[k];
    paramsHash = HashParams[params];

    Print["Case ", k, "..."];
    data = Quiet @ Check[GetExampleData[k], $Failed];
    If[data === $Failed,
      Print["  SKIPPED (GetExampleData failed)"];
      Continue[]
    ];
    If[params =!= {}, data = data /. params];

    eqs = Quiet @ Check[DataToEquations[data], $Failed];
    If[eqs === $Failed,
      Print["  SKIPPED (DataToEquations failed)"];
      AppendTo[entries,
        <|"case" -> ToString[k], "params_hash" -> paramsHash,
          "status" -> "D2E_FAILED", "solution_hash" -> Null, "num_keys" -> 0|>
      ];
      Continue[]
    ];

    ClearSolveCache[];
    result = Quiet @ Check[CriticalCongestionSolver[eqs], $Failed];
    If[result === $Failed,
      Print["  SKIPPED (CriticalCongestionSolver failed)"];
      AppendTo[entries,
        <|"case" -> ToString[k], "params_hash" -> paramsHash,
          "status" -> "SOLVER_FAILED", "solution_hash" -> Null, "num_keys" -> 0|>
      ];
      Continue[]
    ];

    status = Lookup[result, "Status", "Unknown"];
    asso = Lookup[result, "AssoCritical", Null];
    solHash = HashSolution[asso];
    nKeys = If[AssociationQ[asso], Length[asso], 0];

    Print["  Status: ", status, " | Keys: ", nKeys, " | Hash: ",
      If[StringQ[solHash], StringTake[solHash, 12] <> "...", "None"]];

    AppendTo[entries,
      <|"case" -> ToString[k], "params_hash" -> paramsHash,
        "status" -> status, "solution_hash" -> solHash, "num_keys" -> nKeys|>
    ];
  ],
  {k, cases}
];

(* --- Merge with existing entries and write JSON --- *)

outPath = FileNameJoin[{Directory[], "MFGraphs", "Tests", "reference_solutions.json"}];

(* Load existing entries and merge: keep entries for cases not in this run *)
existingEntries = If[FileExistsQ[outPath],
  Quiet @ Check[
    Lookup[Import[outPath, "RawJSON"], "cases", {}],
    {}
  ],
  {}
];
newCaseKeys = ToString /@ cases;
mergedEntries = Join[
  Select[existingEntries, !MemberQ[newCaseKeys, #["case"]] &],
  entries
];

(* Build JSON structure *)
output = <|
  "generated" -> DateString[],
  "commit" -> Quiet @ Check[
    StringTrim @ RunProcess[{"git", "rev-parse", "--short", "HEAD"}, "StandardOutput"],
    "unknown"
  ],
  "mathematica_version" -> ToString[$VersionNumber],
  "cases" -> mergedEntries
|>;

Export[outPath, output, "JSON"];
Print["\nWrote ", Length[mergedEntries], " entries (", Length[entries], " new/updated) to ", outPath];
