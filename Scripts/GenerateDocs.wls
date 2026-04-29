#!/usr/bin/env wolframscript
(* GenerateDocs.wls -- Exports MFGraphs usage metadata to API_REFERENCE.md *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = DirectoryName[DirectoryName[$InputFileName]];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];

Print["Generating API Reference..."];

$PublicContexts = {
  "primitives`",
  "scenarioTools`",
  "examples`",
  "unknownsTools`",
  "systemTools`",
  "solversTools`",
  "orchestrationTools`",
  "graphicsTools`"
};

symbols = DeleteDuplicates @ Flatten[Names[# <> "*"] & /@ $PublicContexts];
publicSymbols = Select[
  symbols,
  !StringContainsQ[#, "Private`"] && !StringEndsQ[#, "$"] &
];

mdContent = {
  "# MFGraphs API Reference",
  "",
  "This document is automatically generated from the Wolfram Language package usage metadata.",
  ""
};

Do[
  symName = Last @ StringSplit[s, "`"];
  usage = Information[Symbol[s], "Usage"];
  If[StringQ[usage] && !StringStartsQ[usage, "False is the symbol for"],
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
