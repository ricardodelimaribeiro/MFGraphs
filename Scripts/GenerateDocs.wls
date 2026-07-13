#!/usr/bin/env wolframscript
(* GenerateDocs.wls -- Exports MFGraphs usage metadata to API_REFERENCE.md *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = DirectoryName[DirectoryName[$InputFileName]];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];

(* Shared usage-arity lint (also enforced by MFGraphs/Tests/usage-arity.mt
   in the fast suite): doc generation refuses to publish a call pattern
   that contradicts its definition's arity. *)
Get[FileNameJoin[{$RepoRoot, "Scripts", "UsageArityLint.wls"}]];

Print["Generating API Reference..."];

report = usageArityReport[];

Scan[Print["Arity lint warning: ", #] &, report["Warnings"]];
If[report["Issues"] =!= {},
  Print["ARITY LINT FAILED -- usage strings disagree with definitions:"];
  Scan[Print["  ", #] &, report["Issues"]];
  Print["API_REFERENCE.md was NOT regenerated. Fix the ::usage strings (or definitions) first."];
  Exit[1]
];
If[report["Checked"] === 0,
  Print["ARITY LINT FAILED -- zero symbols were checked; the definition introspection no longer matches the package's definition style."];
  Print["API_REFERENCE.md was NOT regenerated."];
  Exit[1]
];
Print["Arity lint: ", report["Checked"], " symbols checked (",
  report["NoPattern"], " without call patterns, ",
  report["NoDefs"], " without introspectable definitions",
  If[report["Ellipsis"] > 0,
    ", " <> ToString[report["Ellipsis"]] <> " ellipsis shorthand group" <>
      If[report["Ellipsis"] === 1, "", "s"] <> " skipped", ""],
  "). All call patterns match."];

mdContent = {
  "# MFGraphs API Reference",
  "",
  "This document is automatically generated from the Wolfram Language package usage metadata.",
  ""
};

Do[
  If[StringQ[entry[[3]]],
    AppendTo[mdContent, "## " <> entry[[2]]];
    AppendTo[mdContent, ""];
    AppendTo[mdContent, entry[[3]]];
    AppendTo[mdContent, ""];
  ],
  {entry, report["Entries"]}
];

exportPath = FileNameJoin[{$RepoRoot, "API_REFERENCE.md"}];
Export[exportPath, StringRiffle[mdContent, "\n"], "Text"];

Print["Documentation exported to: ", exportPath];
(* :!CodeAnalysis::EndBlock:: *)
