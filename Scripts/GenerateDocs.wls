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
  "utilities`",
  "scenarioTools`",
  "examples`",
  "unknownsTools`",
  "systemTools`",
  "solversTools`",
  "orchestrationTools`",
  "graphicsTools`"
};

(* `Names[ctx<>"*"]` already matches `$`-prefixed symbols, so a single pass
   suffices. ($-symbols are surfaced correctly by `readUsage` below, not here.) *)
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

(* Symbols with values (e.g. $MFGraphsVerbose = False) cause
   `Information[Symbol[s], "Usage"]` to read the value's usage instead of the
   symbol's. Use `MessageName` on a held symbol so the usage is read off the
   symbol metadata directly. *)
readUsage[name_String] := ToExpression[name <> "::usage", InputForm];

Do[
  symName = Last @ StringSplit[s, "`"];
  usage = readUsage[s];
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
