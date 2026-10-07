Get[FileNameJoin[{$RepoRoot, "Scripts", "ExactArtifactIntegrity.wls"}]];

Test[
    artifactContentReport[<|"Held" -> HoldComplete[1 + 1]|>,
        <|"Held" -> HoldComplete[2]|>]["Integrity"],
    "ContentMismatch", TestID -> "artifact: held expressions preserve unevaluated structure"
]
Test[
    Module[{counter = 0, payload, report},
        payload = <|"Held" -> HoldComplete[counter++]|>;
        report = artifactContentReport[payload, payload];
        {report["Integrity"], counter}],
    {"ExactContentMatch", 0}, TestID -> "artifact: held metadata has no comparison side effects"
]
Test[
    Block[{heldParameter = 3},
        artifactContentReport[<|"Function" -> Function[heldParameter, heldParameter + 1]|>,
            <|"Function" -> Function[heldParameter, 4]|>]["Integrity"]],
    "ContentMismatch", TestID -> "artifact: local function parameter does not evaluate during comparison"
]

Test[
    With[{payload = <|"Hamiltonian" -> <|"G" -> Function[z, -1/z]|>|>},
        artifactContentReport[payload, BinaryDeserialize[BinarySerialize[payload]]]["Integrity"]],
    "ExactContentMatch", TestID -> "artifact: held function metadata compares without evaluation messages"
]
Test[
    artifactContentReport[<|"G" -> Function[z, -1/z]|>,
        <|"G" -> Function[z, -2/z]|>]["Integrity"],
    "ContentMismatch", TestID -> "artifact: changed Hamiltonian function is detected"
]

Test[
    Module[{payload, reloaded},
        payload = <|"Graph" -> Graph[{1 <-> 2}], "Solution" -> {j[1, 2] -> 1/3},
            "Assumptions" -> <|"Alpha" -> 1, "Entries" -> {{1, 1/3}}, "Domain" -> Reals|>|>;
        reloaded = BinaryDeserialize[BinarySerialize[payload]];
        artifactContentReport[payload, reloaded]["Integrity"]
    ], "ExactContentMatch", TestID -> "artifact: graph and all mathematical metadata survive WXF"
]
Test[
    With[{original = <|"Solution" -> {j[1, 2] -> 1/3}, "Domain" -> Reals|>,
          changed = <|"Solution" -> {j[1, 2] -> 1/3 + 1/10^30}, "Domain" -> Reals|>},
        artifactContentReport[original, changed]["Integrity"]
    ], "ContentMismatch", TestID -> "artifact: detects arbitrarily small exact solution change"
]
Test[
    artifactContentReport[<|"Residual" -> j[1, 2, 3] >= 0|>,
        <|"Residual" -> j[1, 2, 3] > 0|>]["Integrity"],
    "ContentMismatch", TestID -> "artifact: detects lost boundary of parametric domain"
]
Test[
    artifactContentReport[<|"Alpha" -> 1, "Domain" -> Reals|>,
        <|"Alpha" -> 1/2, "Domain" -> Reals|>]["Integrity"],
    "ContentMismatch", TestID -> "artifact: detects changed model assumption"
]
Test[
    artifactContentReport[<|"Solution" -> {primitives`j[1, 2] -> 1}|>,
        <|"Solution" -> {otherContext`j[1, 2] -> 1}|>]["Integrity"],
    "ContentMismatch", TestID -> "artifact: preserves symbol contexts"
]
Test[
    artifactContentReport[<|"Graph" -> Graph[{1 <-> 2}, EdgeWeight -> {1}]|>,
        <|"Graph" -> Graph[{1 <-> 2}, EdgeWeight -> {2}]|>]["Integrity"],
    "ContentMismatch", TestID -> "artifact: detects changed edge metadata"
]
Test[
    artifactContentReport[<|"Adjacency" -> SparseArray[{{1, 2} -> 1}, {2, 2}]|>,
        <|"Adjacency" -> SparseArray[{{1, 2} -> 1}, {3, 3}]|>]["Integrity"],
    "ContentMismatch", TestID -> "artifact: preserves sparse dimensions"
]
Test[
    artifactContentReport[<|"Solution" -> {}, "Status" -> "Timeout"|>,
        <|"Solution" -> {}|>]["Integrity"],
    "ContentMismatch", TestID -> "artifact: missing metadata cannot pass"
]
