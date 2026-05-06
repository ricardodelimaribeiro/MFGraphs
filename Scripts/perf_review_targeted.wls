#!/usr/bin/env wolframscript
(* Usage: wolframscript -file Scripts/perf_review_targeted.wls *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

SetAttributes[measure, HoldAll];

measure[label_String, expr_] := Module[{mem0, max0, t, result, mem1, max1},
  mem0 = MemoryInUse[];
  max0 = MaxMemoryUsed[];
  {t, result} = AbsoluteTiming[expr];
  mem1 = MemoryInUse[];
  max1 = MaxMemoryUsed[];
  <|
    "Label" -> label,
    "TimeSeconds" -> N[t],
    "NetMemoryBytes" -> (mem1 - mem0),
    "PeakExtraBytes" -> (max1 - max0),
    "Result" -> result
  |>
];

safeReduce[sys_] := TimeConstrained[reduceSystem[sys], 20, $Aborted];

scenarioSpecs = {
  <|"Name" -> "chain_3", "Build" -> Function[{}, gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}]], "RunReduce" -> True|>,
  <|"Name" -> "chain3_mid_sc", "Build" -> Function[{}, gridScenario[{3}, {{2, 80.0}}, {{1, 0.0}, {3, 10.0}}, {{1,2,3,2},{3,2,1,2}}]], "RunReduce" -> True|>,
  <|"Name" -> "attraction_11", "Build" -> Function[{}, getExampleScenario[11, {{1, 100.0}}, {{4, 0.0}}]], "RunReduce" -> True|>,
  <|"Name" -> "jamarat_20", "Build" -> Function[{}, getExampleScenario[20, {{1, 60.0}, {2, 40.0}}, {{7, 0.0}, {8, 1.0}, {9, 2.0}}]], "RunReduce" -> False|>
};

runCase[spec_Association] := Module[
  {name, sM, uM, sysM, redM, sys, result, topology},
  name = spec["Name"];
  sM = measure["scenario", spec["Build"][]];
  If[FailureQ[sM["Result"]], Return[<|"Case" -> name, "ScenarioFailure" -> ToString[sM["Result"], InputForm]|>]];
  uM = measure["symbolicUnknowns", makeSymbolicUnknowns[sM["Result"]]];
  If[FailureQ[uM["Result"]], Return[<|"Case" -> name, "UnknownsFailure" -> ToString[uM["Result"], InputForm]|>]];
  sysM = measure["system", makeSystem[sM["Result"], uM["Result"]]];
  If[FailureQ[sysM["Result"]], Return[<|"Case" -> name, "SystemFailure" -> ToString[sysM["Result"], InputForm]|>]];
  sys = sysM["Result"];
  topology = scenarioData[sM["Result"], "Topology"];
  redM = If[TrueQ[spec["RunReduce"]], measure["reduce", safeReduce[sys]], <|"TimeSeconds" -> "Skipped", "PeakExtraBytes" -> "Skipped", "Result" -> "Skipped"|>];
  result = redM["Result"];
  <|
    "Case" -> name,
    "ScenarioTimeSeconds" -> sM["TimeSeconds"],
    "ScenarioPeakExtraBytes" -> sM["PeakExtraBytes"],
    "UnknownsTimeSeconds" -> uM["TimeSeconds"],
    "UnknownsPeakExtraBytes" -> uM["PeakExtraBytes"],
    "SystemTimeSeconds" -> sysM["TimeSeconds"],
    "SystemPeakExtraBytes" -> sysM["PeakExtraBytes"],
    "ReduceTimeSeconds" -> redM["TimeSeconds"],
    "ReducePeakExtraBytes" -> redM["PeakExtraBytes"],
    "ScenarioByteCount" -> ByteCount[sM["Result"]],
    "UnknownsByteCount" -> ByteCount[uM["Result"]],
    "SystemByteCount" -> ByteCount[sys],
    "HalfPairs" -> Length[Lookup[topology, "HalfPairs", {}]],
    "AuxPairs" -> Length[Lookup[topology, "AuxPairs", {}]],
    "AuxTriples" -> Length[Lookup[topology, "AuxTriples", {}]],
    "Js" -> Length[symbolicUnknownsData[uM["Result"], "Js"]],
    "Jts" -> Length[symbolicUnknownsData[uM["Result"], "Jts"]],
    "ReduceResultKind" -> Which[
      result === $Aborted, "TimedOut",
      result === False, "False",
      ListQ[result], "Rules",
      AssociationQ[result], "Association",
      result === "Skipped", "Skipped",
      True, Head[result] // ToString
    ]
  |>
];

results = runCase /@ scenarioSpecs;
Export[FileNameJoin[{$RepoRoot, "Scripts", "perf_review_targeted_results.json"}], results, "JSON"];
Print[ExportString[results, "JSON", "Compact" -> False]];
