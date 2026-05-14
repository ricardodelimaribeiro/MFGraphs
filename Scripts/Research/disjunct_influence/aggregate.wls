#!/usr/bin/env wolframscript
(* PR 2: aggregate per-scenario JSON from measure_branching.wls and
   measure_structural_role.wls into a single Markdown summary.

   Two sections:
     - Branching (substituted-form granularity): top-N most-branching
       disjuncts per scenario; cross-scenario distribution of pruneRate.
     - Structural (abstract-complementarity granularity): per-feature
       conditional probabilities P(outcome | feature bin), looking for
       a topology rule that fits all cases.

   Reads from:
     Scripts/Research/results/branching/<key>.json
     Scripts/Research/results/structural/<key>.json
   Writes to:
     Scripts/Research/results/findings.md *)

$RepoRoot = ParentDirectory[ParentDirectory[ParentDirectory[DirectoryName[$InputFileName]]]];
$resultsDir = FileNameJoin[{$RepoRoot, "Scripts", "Research", "results"}];
$branchingDir  = FileNameJoin[{$resultsDir, "branching"}];
$structuralDir = FileNameJoin[{$resultsDir, "structural"}];
$outFile = FileNameJoin[{$resultsDir, "findings.md"}];

(* ---- Loaders ---- *)

loadAll[dir_] :=
    Module[{files},
        files = FileNames["*.json", dir];
        Association @ Map[
            FileBaseName[#] -> Quiet[Import[#, "RawJSON"]] &,
            files
        ]
    ];

$branchingData  = If[DirectoryQ[$branchingDir],  loadAll[$branchingDir],  <||>];
$structuralData = If[DirectoryQ[$structuralDir], loadAll[$structuralDir], <||>];

Print["Loaded branching: ", Length[$branchingData], " scenarios"];
Print["Loaded structural: ", Length[$structuralData], " scenarios"];

(* ---- Branching aggregate ---- *)

branchingTopN[key_, n_:5] :=
    Module[{rows},
        rows = Lookup[$branchingData[key], "perDisjunct", {}];
        Take[rows, UpTo[n]]
    ];

branchingSummaryRow[key_] :=
    Module[{data, totals, conservation, top},
        data = $branchingData[key];
        totals = Lookup[data, "totals", <||>];
        conservation = Lookup[data, "conservation", <||>];
        top = branchingTopN[key, 1];
        <|
            "Scenario" -> key,
            "Status"   -> Lookup[data, "status", "?"],
            "Started"  -> Lookup[totals, "BranchesStarted", "?"],
            "Unique"   -> Lookup[totals, "UniqueDisjuncts", "?"],
            "MaxDepth" -> Lookup[totals, "MaxRecursionDepth", "?"],
            "TopShare" -> If[Length[top] > 0,
                NumberForm[N[top[[1, "started"]] / Max[1, Lookup[totals, "BranchesStarted", 1]]], {4, 3}],
                "0"]
        |>
    ];

(* ---- Structural aggregate ---- *)

(* Bin thresholds derived from the actual value ranges observed across
   all 38 scenarios (computed in PR 2 refinement):
     edgeBetweenness: min 2, q1 73, median 285, q3 740, max 2770
     boundaryDistance: 0..14 (grids reach 14 hops)
     edgeCost: discrete small integers 0..5+
     junctionDegree: 1..6
   The aggregate uses these bins; PR 4 may iterate further. *)
binFeature["junctionDegree", v_?NumericQ] :=
    Which[v <= 2, "small(<=2)", v <= 4, "med(3-4)", True, "large(>=5)"];
binFeature["boundaryDistance", v_?NumericQ] :=
    Which[v == 0, "0", v == 1, "1", v == 2, "2", v == 3, "3",
          v <= 5, "4-5", v <= 9, "6-9", True, "10+"];
binFeature["edgeBetweenness", v_?NumericQ] :=
    Which[v < 50, "low(<50)", v < 200, "mid(50-200)",
          v < 800, "high(200-800)", True, "very-high(>=800)"];
binFeature["edgeCost", v_?NumericQ] :=
    Which[v == 0, "0", v <= 1, "(0,1]", v <= 5, "(1,5]", True, "(5,inf)"];
binFeature[_, _] := "missing";

structuralRowsAcrossScenarios[cachedOnly_:False] :=
    Flatten[
        Function[k,
            Module[{data = $structuralData[k]},
                If[cachedOnly && !TrueQ[Lookup[data, "hasCachedSolution", False]],
                    {},
                    With[{rows = Lookup[data, "rows", {}]},
                        <|"scenario" -> k, "row" -> #|> & /@ rows
                    ]
                ]
            ]
        ] /@ Keys[$structuralData],
        1
    ];

featureOutcomeTable[allRows_, feature_] :=
    Module[{groups, counts, table},
        groups = GroupBy[allRows, binFeature[feature, #["row"][feature]] &];
        table = Function[bin,
            Module[{outcomes = #["row"]["outcome"] & /@ groups[bin]},
                <|
                    "feature"  -> feature,
                    "bin"      -> bin,
                    "n"        -> Length[outcomes],
                    "Inactive" -> Count[outcomes, "Inactive"],
                    "Active"   -> Count[outcomes, "Active"],
                    "Ambiguous"-> Count[outcomes, "Ambiguous"],
                    "P_Inactive" -> NumberForm[
                        N[Count[outcomes, "Inactive"] / Max[1, Length[outcomes]]],
                        {4, 3}],
                    "P_Active"   -> NumberForm[
                        N[Count[outcomes, "Active"] / Max[1, Length[outcomes]]],
                        {4, 3}]
                |>
            ]
        ] /@ Keys[groups];
        SortBy[table, #["bin"] &]
    ];

(* ---- Markdown emission ---- *)

mdRow[assoc_, headers_] := "| " <> StringRiffle[
    Function[h, ToString[Lookup[assoc, h, ""]]] /@ headers,
    " | "] <> " |";

mdTable[rows_List, headers_List] :=
    StringRiffle[
        Join[
            {"| " <> StringRiffle[headers, " | "] <> " |"},
            {"| " <> StringRiffle[ConstantArray["---", Length[headers]], " | "] <> " |"},
            mdRow[#, headers] & /@ rows
        ],
        "\n"
    ];

content = StringJoin[
    "# Disjunct-influence study — aggregate findings\n\n",
    "Generated by `Scripts/Research/disjunct_influence/aggregate.wls`.\n\n",
    "## Branching summary (substituted-form granularity)\n\n",
    "Per scenario: total branches started, unique substituted-form disjunct\n",
    "IDs, max recursion depth, and the share of branching attributable to\n",
    "the single most-branching substituted disjunct (`TopShare`).\n\n",
    If[$branchingData === <||>, "_No branching data; run `measure_branching.wls`._\n",
        mdTable[
            SortBy[branchingSummaryRow /@ Keys[$branchingData], -ToExpression[ToString[#["Started"]]] &],
            {"Scenario", "Status", "Started", "Unique", "MaxDepth", "TopShare"}
        ]
    ],
    "\n\n## Structural-role distribution — cached-only (abstract-complementarity granularity)\n\n",
    "Filtered to scenarios with a cached symbolic solution (so Inactive/Active\n",
    "labels are real, not defaulted). `P_Inactive` close to 1 means the feature\n",
    "is a strong topology-only predictor that the complementarity is inactive\n",
    "at equilibrium; `P_Active` close to 1 is the symmetric finding.\n\n",
    If[$structuralData === <||>, "_No structural data; run `measure_structural_role.wls`._\n",
        Module[{cachedRows = structuralRowsAcrossScenarios[True]},
            StringJoin @@ (
                Function[feature,
                    "### " <> feature <> "\n\n" <>
                    mdTable[
                        featureOutcomeTable[cachedRows, feature],
                        {"bin", "n", "Inactive", "Active", "Ambiguous", "P_Inactive", "P_Active"}
                    ] <> "\n\n"
                ] /@ {"junctionDegree", "boundaryDistance", "edgeBetweenness", "edgeCost"}
            )
        ]
    ],
    "\n## Structural-role distribution — all scenarios (incl. uncached)\n\n",
    "Same features over the full dataset. The Ambiguous bucket is inflated\n",
    "by 5 grids without cached solutions (Grid0505 / Grid0707 / Grid0710 /\n",
    "Grid1010 / Grid1020 — 2,405 disjuncts). Use the cached-only section\n",
    "above for actual conditional probabilities.\n\n",
    If[$structuralData === <||>, "",
        Module[{allRows = structuralRowsAcrossScenarios[False]},
            StringJoin @@ (
                Function[feature,
                    "### " <> feature <> "\n\n" <>
                    mdTable[
                        featureOutcomeTable[allRows, feature],
                        {"bin", "n", "Inactive", "Active", "Ambiguous", "P_Inactive", "P_Active"}
                    ] <> "\n\n"
                ] /@ {"junctionDegree", "boundaryDistance", "edgeBetweenness", "edgeCost"}
            )
        ]
    ],
    "\n## Decision\n\n",
    "_Hand-edit this section after reading the tables above. Three forks:_\n",
    "1. Strong topology rule visible (e.g., `P(Inactive | boundaryDistance >= 3) > 0.95`):\n",
    "   ship a `presolveByTopology` helper.\n",
    "2. Most branching concentrated in one substituted disjunct per scenario:\n",
    "   ship a per-scenario disjunct-priority hint.\n",
    "3. Neither: document the negative result and consider Phase B (block ordering)\n",
    "   as the next test.\n"
];

Export[$outFile, content, "Text"];
Print["Wrote ", $outFile];
