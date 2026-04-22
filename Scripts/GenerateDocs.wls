#!/usr/bin/env wolframscript
(* GenerateDocs.wls -- Exports MFGraphs usage metadata to API_REFERENCE.md *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = "/Users/ribeirrd/Documents/GitHub/MFGraphs";
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];

Print["Generating API Reference..."];

symbols = Names["MFGraphs`*"];
publicSymbols = Select[symbols, (!StringContainsQ[#, "Private`"]) &];

mdContent = {
  "# MFGraphs API Reference",
  "",
  "This document is automatically generated from the Wolfram Language package usage metadata.",
  ""
};

Do[
  symName = StringReplace[s, "MFGraphs`" -> ""];
  usage = Information[Symbol[s], "Usage"];
  If[StringQ[usage],
    AppendTo[mdContent, "## " <> symName];
    AppendTo[mdContent, ""];
    AppendTo[mdContent, usage];
    AppendTo[mdContent, ""];
  ],
  {s, publicSymbols}
];

exportPath = FileNameJoin[{$RepoRoot, "API_REFERENCE.md"}];
Export[exportPath, StringRiffle[mdContent, "\n"], "Text"];

Print["Documentation exported to: ", exportPath];
(* :!CodeAnalysis::EndBlock:: *)
