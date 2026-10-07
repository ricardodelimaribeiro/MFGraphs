#!/usr/bin/env wolframscript
(* Verify each saved file against its recorded checksum BEFORE importing it,
   then recheck soundness and a fresh content-preserving round trip separately.
   Checksums establish correspondence with the recorded file, not mathematical
   validity or authenticity of a separately altered checksum manifest.
   Recheck saved solutions without solving; export readable exact mathematics.
   Usage: wolframscript -file Scripts/AuditCuratedResults.wls selection.json NEW-output-directory
   selection.json maps case names to repo-relative result.wxf paths. *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
$RepoRoot = ParentDirectory[DirectoryName[$InputFileName]];
PrependTo[$Path, $RepoRoot]; Needs["MFGraphs`"];
Get[FileNameJoin[{$RepoRoot, "Scripts", "ExactArtifactIntegrity.wls"}]];
args = Select[Rest[$ScriptCommandLine], !StringEndsQ[#, ".wls"] &];
If[Length[args] =!= 2, Print["Expected selection.json and a NEW output directory"]; Exit[2]];
{selectionPath, outDir} = args;
If[DirectoryQ[outDir], Print["Output directory already exists; refusing overwrite"]; Exit[2]];
CreateDirectory[outDir, CreateIntermediateDirectories -> True];
selection = Import[selectionPath, "RawJSON"];
If[!AssociationQ[selection] || Length[selection] == 0 ||
    !AllTrue[Normal[selection], MatchQ[#, _String -> _String] &],
    Print["Selection must be a nonempty association of case names to artifact paths"]; Exit[2]];
rows = {}; allValid = True;
KeyValueMap[Function[{name, relativePath}, Module[
    {path, payload, saved, fresh, content, readable, roundtrip, row,
     markerPath, marker, expectedHash, actualHash, recordedIntegrity},
    path = FileNameJoin[{$RepoRoot, relativePath}];
    markerPath = FileNameJoin[{DirectoryName[path], "completed.json"}];
    marker = If[FileExistsQ[markerPath], Quiet[Check[Import[markerPath, "RawJSON"], $Failed]], $Failed];
    expectedHash = If[AssociationQ[marker], Lookup[marker, "ArtifactSHA256", Missing["NotRecorded"]],
        Missing["NotRecorded"]];
    actualHash = If[FileExistsQ[path], FileHash[path, "SHA256", "HexString"], Missing["FileAbsent"]];
    recordedIntegrity = Which[
        MissingQ[actualHash], "ArtifactMissing",
        !StringQ[expectedHash], "RecordedChecksumMissing",
        !TrueQ[Lookup[marker, "Complete", False]], "RecordedCompletionMissing",
        expectedHash =!= actualHash, "RecordedChecksumMismatch",
        True, "RecordedChecksumMatch"];
    If[recordedIntegrity =!= "RecordedChecksumMatch",
        row = <|"Case" -> name, "Artifact" -> relativePath,
            "RecordedFileIntegrity" -> recordedIntegrity, "RecordedArtifactSHA256" -> expectedHash,
            "ArtifactSHA256" -> actualHash, "ArtifactIntegrity" -> "NotChecked",
            "RecheckSoundness" -> "NotChecked"|>;
        Export[FileNameJoin[{outDir, name <> ".wl"}], row, "WL"];
        (* Missing objects are rendered explicitly for the JSON audit ledger. *)
        AppendTo[rows, row /. Missing[reason_] :> reason];
        allValid = False;
        Print[name, ": ", recordedIntegrity, "; artifact was not imported"];
        Return[Null, Module]];
    payload = Import[path, "WXF"]; saved = payload["Validation"];
    fresh = exactSolutionReport[payload["System"], payload["Solution"],
        "TimeLimit" -> 10, "CheckCompleteness" -> False];
    (* This comparison establishes content preservation of this fresh round trip;
       the recorded checksum comparison above establishes historical continuity. *)
    roundtrip = BinaryDeserialize[BinarySerialize[payload]];
    content = artifactContentReport[payload, roundtrip];
    readable = <|"Artifact" -> relativePath, "SHA256" -> actualHash,
        "RecordedArtifactSHA256" -> expectedHash, "RecordedFileIntegrity" -> recordedIntegrity,
        "Specification" -> payload["Specification"], "Solution" -> payload["Solution"],
        "ParameterDomain" -> saved["ParameterDomain"], "FreeVariables" -> saved["FreeVariables"],
        "OriginalConstraints" -> saved["OriginalConstraints"],
        "SavedValidation" -> KeyDrop[saved, {"OriginalConstraints", "SolutionExpression"}],
        "FreshSoundnessRecheck" -> KeyDrop[fresh, {"OriginalConstraints", "SolutionExpression"}],
        "ContentIntegrity" -> content|>;
    Export[FileNameJoin[{outDir, name <> ".wl"}], readable, "WL"];
    row = Join[payload["Summary"], <|"Artifact" -> relativePath,
        "ArtifactSHA256" -> readable["SHA256"],
        "RecordedFileIntegrity" -> recordedIntegrity,
        "ArtifactIntegrity" -> content["Integrity"],
        "SoundnessSeconds" -> saved["SoundnessSeconds"],
        "CompletenessSeconds" -> saved["CompletenessSeconds"],
        "FreeVariableCount" -> Length[saved["FreeVariables"]],
        "RecheckSoundness" -> fresh["Soundness"],
        "RecheckSeconds" -> fresh["SoundnessSeconds"]|>];
    AppendTo[rows, row];
    allValid = allValid && saved["Soundness"] === "Proved" && fresh["Soundness"] === "Proved" &&
        recordedIntegrity === "RecordedChecksumMatch" &&
        content["Integrity"] === "ExactContentMatch";
    Print[name, ": saved=", saved["Soundness"], " recheck=", fresh["Soundness"],
        " completeness=", saved["Completeness"], " integrity=", content["Integrity"]];
]], selection];
Export[FileNameJoin[{outDir, "results.json"}], rows, "RawJSON"];
Export[FileNameJoin[{outDir, "completed.json"}], <|"Complete" -> allValid,
    "ExpectedCases" -> Length[selection], "ActualCases" -> Length[rows],
    "WolframVersion" -> $Version,
    "ValidatorSHA256" -> FileHash[FileNameJoin[{$RepoRoot, "MFGraphs", "validationTools.wl"}], "SHA256", "HexString"],
    "IntegrityScriptSHA256" -> FileHash[FileNameJoin[{$RepoRoot, "Scripts", "ExactArtifactIntegrity.wls"}], "SHA256", "HexString"],
    "AuditScriptSHA256" -> FileHash[$InputFileName, "SHA256", "HexString"]|>, "RawJSON"];
Exit[If[allValid, 0, 1]];
