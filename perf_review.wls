Needs["MFGraphs`"];

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

scenarioSpecs = {
  <|"Name" -> "chain_3", "Build" -> Function[{}, gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}]]|>,
  <|"Name" -> "y_case_8", "Build" -> Function[{}, getExampleScenario[8, {{1, 80.0}}, {{3, 0.0}, {4, 10.0}}]]|>,
  <|"Name" -> "attraction_11", "Build" -> Function[{}, getExampleScenario[11, {{1, 100.0}}, {{4, 0.0}}]]|>,
  <|"Name" -> "jamarat_20", "Build" -> Function[{}, getExampleScenario[20, {{1, 60.0}, {2, 40.0}}, {{7, 0.0}, {8, 1.0}, {9, 2.0}}]]|>
};

runCase[spec_Association] := Module[
  {name, sM, uM, sysM, redM, sys, result, out, topology},
  name = spec["Name"];
  sM = measure["scenario", spec["Build"][]];
  If[FailureQ[sM["Result"]],
    Return[<|"Case" -> name, "ScenarioFailure" -> ToString[sM["Result"], InputForm]|>]
  ];
  uM = measure["unknowns", makeUnknowns[sM["Result"]]];
  If[FailureQ[uM["Result"]],
    Return[<|"Case" -> name, "UnknownsFailure" -> ToString[uM["Result"], InputForm]|>]
  ];
  sysM = measure["system", makeSystem[sM["Result"], uM["Result"]]];
  If[FailureQ[sysM["Result"]],
    Return[<|"Case" -> name, "SystemFailure" -> ToString[sysM["Result"], InputForm]|>]
  ];
  sys = sysM["Result"];
  topology = scenarioData[sM["Result"], "Topology"];
  redM = measure["reduce", reduceSystem[sys]];
  result = redM["Result"];
  out = <|
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
    "AuxPairs" -> Length[Lookup[topology, "AuxPairs", {}]],
    "AuxTriples" -> Length[Lookup[topology, "AuxTriples", {}]],
    "Js" -> Length[unknownsData[uM["Result"], "Js"]],
    "Jts" -> Length[unknownsData[uM["Result"], "Jts"]],
    "ReduceResultKind" -> Which[
      result === False, "False",
      ListQ[result], "Rules",
      AssociationQ[result], "Association",
      True, Head[result] // ToString
    ]
  |>; 
  out
];

results = runCase /@ scenarioSpecs;
Export["perf_review_results.json", results, "JSON"];
Print[ExportString[results, "JSON", "Compact" -> False]];
