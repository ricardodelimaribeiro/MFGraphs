#!/usr/bin/env wolframscript
(* Read-only diagnosis of original vs WXF-reloaded mathematical content.
   Creates a new evidence directory; never overwrites the failing run. *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
$RepoRoot = ParentDirectory[DirectoryName[$InputFileName]];
PrependTo[$Path, $RepoRoot]; Needs["MFGraphs`"];
Get[FileNameJoin[{$RepoRoot, "MFGraphs", "validationTools.wl"}]];
outDir = FileNameJoin[{$RepoRoot, "Results", "exact-foundation",
    DateString[{"Year", "Month", "Day", "T", "Hour", "Minute", "Second"}] <>
        "-roundtrip-diagnosis-" <> StringTake[CreateUUID[], 8]}];
CreateDirectory[outDir];
s = gridScenario[{2}, {{1, 10}}, {{2, 0}}];
sys = makeSystem[s]; sol = dnfReduceSystem[sys];
validation = exactSolutionReport[sys, sol, "TimeLimit" -> 2];
original = <|"Scenario" -> s, "System" -> sys, "Solution" -> sol,
    "Validation" -> validation, "Assumptions" -> <|"Domain" -> Reals, "Alpha" -> 1|>|>;
path = FileNameJoin[{outDir, "reproduced-failure.wxf"}];
Export[path, original, "WXF"];
reloaded = Import[path, "WXF"];

graphView[g_Graph] := <|"Vertices" -> VertexList[g], "Edges" -> EdgeList[g],
    "Options" -> Options[g], "InputForm" -> ToString[g, InputForm]|>;
diffs[a_, b_, path_List : {}] := Module[{children}, Which[
    SameQ[a, b], {},
    Head[a] === Graph && Head[b] === Graph,
        {<|"Path" -> path, "Kind" -> "GraphObject", "AtomQ" -> {AtomQ[a], AtomQ[b]},
          "VertexListSame" -> SameQ[VertexList[a], VertexList[b]],
          "EdgeListSame" -> SameQ[EdgeList[a], EdgeList[b]],
          "OptionsSame" -> SameQ[Options[a], Options[b]],
          "InputFormSame" -> SameQ[ToString[a, InputForm], ToString[b, InputForm]],
          "Original" -> graphView[a], "Reloaded" -> graphView[b]|>},
    AssociationQ[a] && AssociationQ[b] && Keys[a] === Keys[b],
        children = Flatten[KeyValueMap[diffs[#2, b[#1], Append[path, #1]] &, a], 1];
        If[children === {}, {<|"Path" -> path, "Kind" -> "AssociationRepresentation",
            "NormalSameQ" -> SameQ[Normal[a], Normal[b]],
            "InputFormSame" -> SameQ[ToString[a, InputForm], ToString[b, InputForm]]|>}, children],
    Head[a] === Head[b] && !AtomQ[a] && !AtomQ[b] && Length[a] === Length[b],
        children = Flatten[MapThread[diffs[#1, #2, Append[path, #3]] &,
            {List @@ a, List @@ b, Range[Length[a]]}], 1];
        If[children === {}, {<|"Path" -> path, "Kind" -> "ExpressionRepresentation",
            "Head" -> Head[a], "InputFormSame" -> SameQ[ToString[a, InputForm], ToString[b, InputForm]]|>}, children],
    True, {<|"Path" -> path, "Kind" -> "ExpressionDifference", "Original" -> a, "Reloaded" -> b|>}
]];
differences = diffs[original, reloaded];
canonicalGraphs[expr_] := expr /. g_Graph :> graphView[g];
evidence = <|
    "WolframVersion" -> $Version,
    "RawSameQ" -> SameQ[original, reloaded],
    "FieldSameQ" -> AssociationMap[SameQ[original[#], reloaded[#]] &, Keys[original]],
    "Differences" -> differences,
    "MinimalGraphAssociation" -> With[{a = <|"Graph" -> Graph[{1 <-> 2}]|>},
        With[{b = BinaryDeserialize[BinarySerialize[a]]},
            <|"AssociationSameQ" -> SameQ[a, b], "NormalSameQ" -> SameQ[Normal[a], Normal[b]],
              "GraphSameQ" -> SameQ[a["Graph"], b["Graph"]]|>]],
    "MinimalScalarAssociation" -> With[{a = <|"Flow" -> 1/3|>},
        SameQ[a, BinaryDeserialize[BinarySerialize[a]]]],
    "GraphContentAndAllOtherFieldsSameQ" -> SameQ[canonicalGraphs[original], canonicalGraphs[reloaded]],
    "OriginalValidation" -> validation,
    "ReloadedValidation" -> exactSolutionReport[reloaded["System"], reloaded["Solution"], "TimeLimit" -> 2]
|>;
Export[FileNameJoin[{outDir, "evidence.wl"}], evidence, "WL"];
Export[FileNameJoin[{outDir, "original-inputform.txt"}], ToString[original, InputForm], "Text"];
Export[FileNameJoin[{outDir, "reloaded-inputform.txt"}], ToString[reloaded, InputForm], "Text"];
Print["EVIDENCE_DIRECTORY=", outDir];
Print["FIELD_SAMEQ=", evidence["FieldSameQ"]];
Print["DIFFERENCES=", InputForm[differences]];
Print["ALL_EXPOSED_CONTENT_SAMEQ=", evidence["GraphContentAndAllOtherFieldsSameQ"]];
Print["MINIMAL_GRAPH_ASSOCIATION=", evidence["MinimalGraphAssociation"]];
Print["MINIMAL_SCALAR_ASSOCIATION=", evidence["MinimalScalarAssociation"]];
