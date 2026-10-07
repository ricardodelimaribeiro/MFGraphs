#!/usr/bin/env wolframscript
(* Print exact mathematical artifacts without triggering a fresh solve.
   Usage: wolframscript -file Scripts/InspectCuratedResult.wls path/to/result.wxf *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
$RepoRoot = ParentDirectory[DirectoryName[$InputFileName]];
PrependTo[$Path, $RepoRoot]; Needs["MFGraphs`"];
args = Select[Rest[$ScriptCommandLine], !StringEndsQ[#, ".wls"] &];
If[Length[args] =!= 1, Print["Expected one result.wxf path"]; Exit[2]];
payload = Import[First[args], "WXF"];
If[!AssociationQ[payload] || !KeyExistsQ[payload, "Solution"], Exit[3]];
Print["SUMMARY: ", InputForm[payload["Summary"]]];
Print["SPECIFICATION: ", InputForm[payload["Specification"]]];
Print["EXACT SOLUTION: ", InputForm[payload["Solution"]]];
Print["VALIDATION: ", InputForm[KeyDrop[payload["Validation"], {"OriginalConstraints", "SolutionExpression"}]]];
