#!/usr/bin/env wolframscript
(* Worker for run_curated_exact.py. One fresh kernel, one case/method/sample. *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
$RunnerRoot = ParentDirectory[DirectoryName[$InputFileName]];
args = Rest[$ScriptCommandLine];
args = Select[args, !StringEndsQ[#, ".wls"] &];
If[Length[args] =!= 7, Print["Expected config.json output-directory source-root case method solve-limit validation-limit"]; Exit[2]];
{configFile, outDir, sourceRoot, caseName, method, solveLimitText, validationLimitText} = args;
config = Import[configFile, "RawJSON"];
If[!AssociationQ[config] || !AssociationQ[config["execution_sources"]] ||
    !MemberQ[Values[config["execution_sources"]], sourceRoot] || config["runner_root"] =!= $RunnerRoot,
    Print["Worker paths do not match the frozen run manifest"]; Exit[2]];
solveLimit = ToExpression[solveLimitText]; validationLimit = ToExpression[validationLimitText];
PrependTo[$Path, sourceRoot];
PrependTo[$Path, FileNameJoin[{sourceRoot, "MFGraphs"}]];
(* Load the same frozen independent validator for every source variant before
   the selected package can load a historical validator of its own. Explicit
   package loading also avoids silently falling back to an installed paclet. *)
Get[FileNameJoin[{$RunnerRoot, "MFGraphs", "validationTools.wl"}]];
Get[FileNameJoin[{sourceRoot, "MFGraphs", "MFGraphs.wl"}]];
$MFGraphsVerbose = False; $HistoryLength = 0;
$RecursionLimit = 16384; $IterationLimit = 1000000;
Get[FileNameJoin[{$RunnerRoot, "Scripts", "CuratedExactScenarios.wls"}]];
Get[FileNameJoin[{$RunnerRoot, "Scripts", "ExactArtifactIntegrity.wls"}]];

runMethod[sys_, "dnf"] := dnfReduceSystem[sys];
runMethod[sys_, "reduce"] := reduceSystem[sys];
runMethod[sys_, "lex"] := Block[{solversTools`Private`$optimizedDNFVarThreshold = Infinity}, activeSetReduceSystem[sys]];
runMethod[sys_, "block-edge"] := activeSetReduceSystem[sys, "DisjunctOrdering" -> "Block-Edge"];
runMethod[sys_, "linear-net"] := linearNetReduceSystem[sys];

(* Identical generic warm-up in every fresh kernel. Never retrieve solveScenario
   memoized answers. Each DNF call's equality cache is invocation-local. *)
{warmTime, warmResult} = AbsoluteTiming[TimeConstrained[
    dnfReduceSystem[makeSystem[gridScenario[{2}, {{1, 10}}, {{2, 0}}]]], 5, $TimedOut]];
{buildTime, sys} = AbsoluteTiming[TimeConstrained[
    scen = curatedExactScenario[caseName]; makeSystem[scen], 15, $TimedOut]];
If[!mfgSystemQ[sys], Print["BUILD FAILED ", InputForm[sys]]; Exit[3]];
{solveTime, sol} = AbsoluteTiming[CheckAbort[TimeConstrained[
    runMethod[sys, method], solveLimit, $TimedOut], $Aborted]];
Print[caseName, " ", method, " solve seconds: ", solveTime];
{validationTime, validation} = AbsoluteTiming[TimeConstrained[
    validationTools`exactSolutionReport[sys, sol, "TimeLimit" -> validationLimit],
    2 validationLimit + 5, <|"ResultType" -> "Unresolved", "Soundness" -> "Timeout", "Completeness" -> "NotChecked"|>]];
{preprocessTime, preprocess} = AbsoluteTiming[TimeConstrained[
    solversTools`Private`buildSolverInputs[sys], Min[5, solveLimit], $TimedOut]];
data = systemDataFlatten[sys];
summary = <|
    "Case" -> caseName, "Method" -> method, "SourceRoot" -> sourceRoot,
    "WolframVersion" -> $Version, "SystemID" -> $SystemID,
    "ProcessorCount" -> $ProcessorCount, "MachineType" -> $MachineType,
    "WarmupSeconds" -> warmTime, "BuildSeconds" -> buildTime,
    "SolveSeconds" -> solveTime, "ValidationSeconds" -> validationTime,
    "PreprocessDiagnosticSeconds" -> preprocessTime,
    "SolveLimit" -> solveLimit, "ValidationLimit" -> validationLimit,
    "Vertices" -> VertexCount[data["Graph"]], "Edges" -> EdgeCount[data["Graph"]],
    "EdgeVariables" -> Length[data["Js"]], "ValueVariables" -> Length[data["Us"]],
    "TransitionVariables" -> Length[data["Jts"]],
    "ResidualVariablesAfterPreprocessing" -> If[ListQ[preprocess], Length[preprocess[[2]]], -1],
    "OrCountAfterPreprocessing" -> If[ListQ[preprocess], Count[preprocess[[1]], _Or, {0, Infinity}], -1],
    "ResultType" -> Lookup[validation, "ResultType", "Unresolved"],
    "Soundness" -> Lookup[validation, "Soundness", "NotChecked"],
    "Completeness" -> Lookup[validation, "Completeness", "NotChecked"]
|>;
payload = <|"Summary" -> summary, "Specification" -> $CuratedExactScenarios[caseName],
    "Scenario" -> scen, "System" -> sys, "Solution" -> sol, "Validation" -> validation,
    "Preprocessing" -> preprocess|>;
artifactPath = FileNameJoin[{outDir, "result.wxf"}];
Export[artifactPath, payload, "WXF"];
reloaded = Import[artifactPath, "WXF"];
integrity = artifactContentReport[payload, reloaded];
Export[FileNameJoin[{outDir, "integrity.wl"}], integrity, "WL"];
If[integrity["Integrity"] =!= "ExactContentMatch", Print["ARTIFACT CONTENT MISMATCH"]; Exit[4]];
summary = Append[summary, "ArtifactIntegrity" -> integrity["Integrity"]];
Export[FileNameJoin[{outDir, "summary.json"}], summary, "RawJSON"];
Export[FileNameJoin[{outDir, "completed.json"}], <|"Complete" -> True,
    "ArtifactSHA256" -> FileHash[artifactPath, "SHA256", "HexString"]|>, "RawJSON"];
Print["COMPLETE ", caseName, " ", method, " ", summary["ResultType"], " ", summary["Soundness"], " ", summary["Completeness"]];
Exit[0];
(* :!CodeAnalysis::EndBlock:: *)
