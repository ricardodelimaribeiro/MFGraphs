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

(* ---------------------------------------------------------------------- *)
(* Arity lint. A usage string that shows name[a, b] while the definition   *)
(* takes five arguments publishes a wrong signature into API_REFERENCE.md. *)
(* Compare every call pattern in each usage string against the arities the *)
(* symbol's DownValues/SubValues accept, and refuse to export on mismatch. *)
(* Symbols with no DownValues/SubValues (inert constructors like j, u,     *)
(* typed heads, $-variables) are skipped and tallied in the summary.       *)
(*                                                                          *)
(* Known limitations: each usage is scanned only for its own symbol's call *)
(* patterns (signatures of OTHER symbols quoted in the text are not        *)
(* checked); context-qualified occurrences (pkg`name[...]) are skipped;    *)
(* sequence-shaped slots (PatternSequence, Longest, Repeated with a count  *)
(* spec) are scored conservatively -- none occur in the package today.     *)
(* ---------------------------------------------------------------------- *)

(* Word-boundary class for WL identifiers. Mirrors the boundary regex in
   Scripts/AuditPublicAPI.py (tests_for_symbol) -- keep the two in sync. *)
symbolNameCharQ[c_String] :=
  StringMatchQ[c, LetterCharacter | DigitCharacter | "$" | "`"];

(* Scan one bracket group starting at openIdx (the "["). Returns
   {tag, posAfterGroup} where tag is the argument count, "Ellipsis" for a
   group containing a top-level ... (documentation shorthand, not an arity
   claim), or "Unterminated" when the group never closes. Depth-counts
   [] {} () <||> and skips string literals, so only top-level commas count. *)
scanBracketGroup[usage_String, openIdx_Integer] :=
  Module[{len = StringLength[usage], i = openIdx + 1, depth = 1, commas = 0,
          sawContent = False, inQuote = False, dotRun = 0, hasEllipsis = False,
          c, c2},
    While[i <= len && depth > 0,
      c  = StringTake[usage, {i}];
      c2 = If[i < len, StringTake[usage, {i + 1}], ""];
      Which[
        inQuote,                 If[c === "\"", inQuote = False],
        c === "\"",              inQuote = True,
        c === "<" && c2 === "|", depth++; i++,
        c === "|" && c2 === ">", depth--; i++,
        c === "[" || c === "{" || c === "(", depth++,
        c === "]" || c === "}" || c === ")", depth--,
        c === "," && depth === 1, commas++
      ];
      If[depth > 0,
        If[!StringMatchQ[c, WhitespaceCharacter], sawContent = True];
        dotRun = If[c === "." && depth === 1 && !inQuote, dotRun + 1, 0];
        If[dotRun >= 3, hasEllipsis = True]
      ];
      i++
    ];
    Which[
      depth > 0,   {"Unterminated", i},
      hasEllipsis, {"Ellipsis", i},
      True,        {If[sawContent, commas + 1, 0], i}
    ]
  ];

(* All call-pattern occurrences of shortName in the usage text. Returns
   <|"Pairs" -> {{firstArity, secondArityOrMissing}, ...},
     "Ellipsis" -> n, "Unterminated" -> n|>.
   Occurrences preceded by a name character (suffix of a longer identifier,
   or context-qualified) or by a double quote (a signature quoted as example
   text) are not arity claims and are skipped. A chained second group
   (name[...][...]) is captured for the curried SubValues check. *)
scanUsageCallPatterns[shortName_String, usage_String] :=
  Module[{len = StringLength[usage], starts, pairs = {}, nEllipsis = 0,
          nUnterminated = 0, prev},
    starts = First /@ StringPosition[usage, shortName <> "["];
    Do[
      Module[{openIdx = start + StringLength[shortName], tag1, next1,
              tag2 = Missing["NoSecondGroup"]},
        If[start > 1,
          prev = StringTake[usage, {start - 1}];
          If[symbolNameCharQ[prev] || prev === "\"", Continue[]]
        ];
        (* name[[...]] is a Part expression in prose, not a call pattern. *)
        If[openIdx < len && StringTake[usage, {openIdx + 1}] === "[",
          Continue[]];
        {tag1, next1} = scanBracketGroup[usage, openIdx];
        Which[
          tag1 === "Unterminated", nUnterminated++,
          tag1 === "Ellipsis",     nEllipsis++,
          True,
            If[next1 <= len && StringTake[usage, {next1}] === "[",
              Module[{t2, n2},
                {t2, n2} = scanBracketGroup[usage, next1];
                Which[
                  t2 === "Unterminated", nUnterminated++,
                  t2 === "Ellipsis",     nEllipsis++,
                  True,                  tag2 = t2
                ]
              ]
            ];
            AppendTo[pairs, {tag1, tag2}]
        ]
      ],
      {start, starts}
    ];
    <|"Pairs" -> pairs, "Ellipsis" -> nEllipsis, "Unterminated" -> nUnterminated|>
  ];

(* {min, max} argument count contributed by a single held argument slot. *)
slotRange[h_Hold] :=
  Replace[h, {
    Hold[Verbatim[Pattern][_, p_]]     :> slotRange[Hold[p]],
    Hold[Verbatim[PatternTest][p_, _]] :> slotRange[Hold[p]],
    Hold[Verbatim[Condition][p_, _]]   :> slotRange[Hold[p]],
    Hold[_OptionsPattern]              :> {0, Infinity},
    Hold[_Optional]                    :> {0, 1},
    Hold[_BlankNullSequence]           :> {0, Infinity},
    Hold[_BlankSequence]               :> {1, Infinity},
    Hold[_RepeatedNull]                :> {0, Infinity},
    Hold[_Repeated]                    :> {1, Infinity},
    _                                  :> {1, 1}
  }];

defArityRange[Hold[args___]] :=
  Module[{ranges = slotRange /@ (List @@ Map[Hold, Hold[args]])},
    {Total[ranges[[All, 1]]], Total[ranges[[All, 2]]]}
  ];

(* Strip whole-LHS Conditions, repeatedly for stacked f[x_] /; c1 /; c2 --
   an unstripped Condition head would itself be miscounted as a 2-argument
   call by the extraction below. *)
stripLhsCondition = Verbatim[HoldPattern][Verbatim[Condition][l_, _]] :> HoldPattern[l];

(* Arity ranges per definition. "Inner" is what a usage's first bracket
   group must satisfy (DownValues plus the inner application of curried
   SubValues); "Outer" is what a chained second group must satisfy. The
   _Symbol head restriction makes curry depth >= 3 fall into the
   uncheckable pool instead of mis-attributing bracket groups. *)
definitionArityRanges[fullName_String] :=
  Module[{dvs, svs},
    dvs = ToExpression[fullName, InputForm, DownValues][[All, 1]] //. stripLhsCondition;
    svs = ToExpression[fullName, InputForm, SubValues][[All, 1]] //. stripLhsCondition;
    <|
      "Inner" -> Join[
        defArityRange /@ Cases[dvs, Verbatim[HoldPattern][_Symbol[args___]] :> Hold[args]],
        defArityRange /@ Cases[svs, Verbatim[HoldPattern][_Symbol[args___][___]] :> Hold[args]]
      ],
      "Outer" ->
        (defArityRange /@ Cases[svs, Verbatim[HoldPattern][_Symbol[___][args___]] :> Hold[args]])
    |>
  ];

(* One shared per-symbol pass feeds both the lint and the export loop, so
   what is linted is exactly what is published. *)
usageEntries = Map[{#, Last @ StringSplit[#, "`"], readUsage[#]} &, publicSymbols];

arityIssues = {};
lintWarnings = {};
checkedCount = 0; noPatternCount = 0; noDefsCount = 0; ellipsisCount = 0;

Do[
  Module[{fullName, shortName, usage, scan, pairs, ranges, badFirst,
          seconds, badSecond},
    {fullName, shortName, usage} = entry;
    If[!StringQ[usage], Continue[]];
    scan = scanUsageCallPatterns[shortName, usage];
    ellipsisCount += scan["Ellipsis"];
    If[scan["Unterminated"] > 0,
      AppendTo[lintWarnings,
        fullName <> ": unterminated call group in usage text (occurrence not checked)"]];
    pairs = scan["Pairs"];
    If[pairs === {}, noPatternCount++; Continue[]];
    ranges = definitionArityRanges[fullName];
    If[ranges["Inner"] === {}, noDefsCount++; Continue[]];
    checkedCount++;
    badFirst = Select[DeleteDuplicates[pairs[[All, 1]]],
      Function[d, NoneTrue[ranges["Inner"], #[[1]] <= d <= #[[2]] &]]];
    seconds = DeleteDuplicates[Cases[pairs, {_, a2_Integer} :> a2]];
    badSecond = If[ranges["Outer"] === {}, {},
      Select[seconds, Function[d, NoneTrue[ranges["Outer"], #[[1]] <= d <= #[[2]] &]]]];
    If[badFirst =!= {},
      AppendTo[arityIssues,
        fullName <> ": usage shows arity " <> ToString[badFirst] <>
          ", definitions accept " <> ToString[ranges["Inner"]]]];
    If[badSecond =!= {},
      AppendTo[arityIssues,
        fullName <> ": usage shows second-group arity " <> ToString[badSecond] <>
          ", curried definitions accept " <> ToString[ranges["Outer"]]]]
  ],
  {entry, usageEntries}
];

Scan[Print["Arity lint warning: ", #] &, lintWarnings];
If[arityIssues =!= {},
  Print["ARITY LINT FAILED -- usage strings disagree with definitions:"];
  Scan[Print["  ", #] &, arityIssues];
  Print["API_REFERENCE.md was NOT regenerated. Fix the ::usage strings (or definitions) first."];
  Exit[1]
];
If[checkedCount === 0,
  Print["ARITY LINT FAILED -- zero symbols were checked; the definition introspection no longer matches the package's definition style."];
  Print["API_REFERENCE.md was NOT regenerated."];
  Exit[1]
];
Print["Arity lint: ", checkedCount, " symbols checked (",
  noPatternCount, " without call patterns, ",
  noDefsCount, " without introspectable definitions",
  If[ellipsisCount > 0,
    ", " <> ToString[ellipsisCount] <> " ellipsis shorthand group" <>
      If[ellipsisCount === 1, "", "s"] <> " skipped", ""],
  "). All call patterns match."];

Do[
  If[StringQ[entry[[3]]],
    AppendTo[mdContent, "## " <> entry[[2]]];
    AppendTo[mdContent, ""];
    AppendTo[mdContent, entry[[3]]];
    AppendTo[mdContent, ""];
  ],
  {entry, usageEntries}
];

exportPath = FileNameJoin[{$RepoRoot, "API_REFERENCE.md"}];
Export[exportPath, StringRiffle[mdContent, "\n"], "Text"];

Print["Documentation exported to: ", exportPath];
(* :!CodeAnalysis::EndBlock:: *)
