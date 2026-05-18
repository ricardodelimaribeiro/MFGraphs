(* ::Package:: *)

(*Quit[]*)


(* ::Title:: *)
(*MFGraphs disjunct-influence study workbook*)


(* ::Subsection:: *)
(*Overview*)


(* ::Text:: *)
(*Companion notebook for the disjunct-influence study (PRs 1-4). Renders the JSON / CSV produced by the four Research scripts:*)
(*	- Scripts/Research/disjunct_influence/measure_branching.wls (PR 1-2): per-substituted-disjunct branching counts*)
(*	- Scripts/Research/disjunct_influence/measure_structural_role.wls (PR 2): per-abstract-complementarity topology features + outcome label*)
(*	- Scripts/Research/disjunct_influence/measure_sensitivity.wls (PR 3): per-disjunct removeDelta + forceResults*)
(*	- Scripts/Research/block_ordering/benchmark_orderings.wls (PR 3): wall time per (scenario, ordering)*)


(* ::Text:: *)
(*Structure:*)
(*	Section 1 \[LongDash] Cross-scenario summary heatmap (P(outcome | feature bin), cached scenarios only)*)
(*	Section 2 \[LongDash] Block-ordering speedup chart (Lex baseline = 1.0)*)
(*	Section 3 \[LongDash] Per-scenario detail (one block per scenario; evaluate independently)*)
(*	Section 4 \[LongDash] Decision document (mirrored to Scripts/Research/results/findings.md)*)


(* ::Subsection:: *)
(*Initialization*)


mfgDir = If[$InputFileName === "",
    NotebookDirectory[],
    DirectoryName[$InputFileName]
];
If[!StringQ[mfgDir] || mfgDir === "", mfgDir = ExpandFileName["."]];
mfgParentDir = ParentDirectory[mfgDir];
If[!MemberQ[$Path, mfgParentDir], PrependTo[$Path, mfgParentDir]];
Needs["MFGraphs`"];

researchResultsDir = FileNameJoin[{mfgParentDir, "Scripts", "Research", "results"}];

loadJSONDir[subdir_] :=
    Module[{dir = FileNameJoin[{researchResultsDir, subdir}], files},
        If[!DirectoryQ[dir], Return[<||>]];
        files = FileNames["*.json", dir];
        Association @ Map[FileBaseName[#] -> Quiet[Import[#, "RawJSON"]] &, files]
    ];

branchingData  = loadJSONDir["branching"];
structuralData = loadJSONDir["structural"];
sensitivityData = loadJSONDir["sensitivity"];

(* Benchmark CSV. *)
benchmarkRaw = Quiet @ Import[
    FileNameJoin[{researchResultsDir, "block_ordering", "benchmark.csv"}], "CSV"];
benchmarkData = If[ListQ[benchmarkRaw] && Length[benchmarkRaw] > 1,
    Module[{headers = First[benchmarkRaw], rows = Rest[benchmarkRaw]},
        AssociationThread[headers, #] & /@ rows],
    {}];

Print["Loaded: ",
    "branching=", Length[branchingData], " ",
    "structural=", Length[structuralData], " ",
    "sensitivity=", Length[sensitivityData], " ",
    "benchmarkRows=", Length[benchmarkData]];


(* ::Section:: *)
(*Section 1 \[LongDash] Cross-scenario topology-feature heatmap*)


(* ::Text:: *)
(*Per-feature conditional probability P(outcome | feature bin) across the 31 cached scenarios. Outcome \[Element] {Inactive, Active, Ambiguous} comes from the cached symbolic solution (see measure_structural_role.wls). High P(Active) on a bin means complementarities at that topology bin overwhelmingly resolve to a non-zero flow at equilibrium; high P(Inactive) is the symmetric finding (the flow pins to zero).*)


(* ::Text:: *)
(*A "fits all cases" topology rule would be a bin with P(outcome) > 0.95 and n > 100. As of PR 2 refinement, no such bin exists; the strongest signal is P(Active | junctionDegree<=2) = 0.80 with n = 218.*)


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

cachedStructuralRows := Flatten[
    Function[k,
        With[{d = structuralData[k]},
            If[!TrueQ[Lookup[d, "hasCachedSolution", False]], {},
                Lookup[d, "rows", {}]]]] /@ Keys[structuralData],
    1];

featureCondTable[feature_String] :=
    Module[{rows = cachedStructuralRows, groups, table},
        groups = GroupBy[rows, binFeature[feature, #[feature]] &];
        table = Function[bin,
            Module[{outcomes = #["outcome"] & /@ groups[bin], n},
                n = Length[outcomes];
                <|
                    "bin"       -> bin,
                    "n"         -> n,
                    "Inactive"  -> Count[outcomes, "Inactive"],
                    "Active"    -> Count[outcomes, "Active"],
                    "Ambiguous" -> Count[outcomes, "Ambiguous"],
                    "P_Inactive" -> N[Count[outcomes, "Inactive"]/Max[1, n]],
                    "P_Active"   -> N[Count[outcomes, "Active"]/Max[1, n]]
                |>]] /@ Keys[groups];
        Dataset[SortBy[table, #["bin"] &]]
    ];

featureHeatmap[feature_String] :=
    Module[{rows = Normal @ featureCondTable[feature], bins, mat},
        bins = #["bin"] & /@ rows;
        mat = {#["P_Inactive"], #["P_Active"]} & /@ rows;
        ArrayPlot[Transpose[mat],
            FrameTicks -> {{{{1, "P_Inactive"}, {2, "P_Active"}}, None},
                           {Thread[{Range[Length[bins]], bins}], None}},
            ColorFunction -> "TemperatureMap",
            ColorFunctionScaling -> False,
            PlotRange -> {0, 1},
            PlotLabel -> Style[feature, 12, Bold],
            ImageSize -> 360,
            FrameLabel -> {"bin", ""}]
    ];


(* Render all four feature heatmaps in a grid. Cached-only crosstab tables
   are also available via featureCondTable[<feature>]. *)
GraphicsGrid[Partition[
    featureHeatmap /@ {"junctionDegree", "boundaryDistance",
                       "edgeBetweenness", "edgeCost"},
    2], ImageSize -> 720]


(* Cached-only crosstab tables \[LongDash] one per feature. *)
Column[
    Function[f,
        Column[{
            Style[f, 13, Bold],
            featureCondTable[f]
        }, Spacings -> 0.3]
    ] /@ {"junctionDegree", "boundaryDistance", "edgeBetweenness", "edgeCost"},
    Spacings -> 1
]


(* ::Section:: *)
(*Section 2 \[LongDash] Block-ordering speedup chart*)


(* ::Text:: *)
(*Wall-time per (scenario, DisjunctOrdering) from PR 3's benchmark. Bars are normalized to Lexicographic = 1.0 within each scenario; bars > 1.0 mean the Block-* ordering is faster. Timeouts (60s) are clamped, so two timing-out orderings show ratio \[TildeTilde] 1.0.*)


benchmarkScenarios := DeleteDuplicates[#["scenario"] & /@ benchmarkData];

scenarioWallTimes[scen_String] :=
    Module[{rows = Select[benchmarkData, #["scenario"] === scen &]},
        Association @ Map[#["ordering"] -> ToExpression[#["wallSec"]] &, rows]
    ];

speedupRow[scen_String] :=
    Module[{w = scenarioWallTimes[scen], lex, base},
        lex = Lookup[w, "Lexicographic", Missing[]];
        If[!NumericQ[lex] || lex <= 0, Return[Missing[]]];
        <|"scenario" -> scen,
          "Lex_s"    -> N[lex],
          "Vertex"   -> N[lex / Lookup[w, "Block-Vertex", lex]],
          "Edge"     -> N[lex / Lookup[w, "Block-Edge",   lex]],
          "SCC"      -> N[lex / Lookup[w, "Block-SCC",    lex]]|>
    ];

speedupRows := DeleteCases[speedupRow /@ benchmarkScenarios, _Missing];

(* Sorted descending by best Block-* speedup. *)
speedupTable := Dataset @ SortBy[speedupRows,
    -Max[#["Vertex"], #["Edge"], #["SCC"]] &];

speedupTable


(* Bar chart \[LongDash] each scenario gets three bars (Vertex / Edge / SCC ratio
   vs Lexicographic). Sorted by best-of-three so the wins cluster at the top. *)
speedupBarChart :=
    Module[{rows = SortBy[Normal[speedupTable],
                          -Max[#["Vertex"], #["Edge"], #["SCC"]] &],
            data, labels},
        data = {#["Vertex"], #["Edge"], #["SCC"]} & /@ rows;
        labels = #["scenario"] & /@ rows;
        BarChart[data,
            ChartLayout -> "Grouped",
            ChartLegends -> {"Block-Vertex", "Block-Edge", "Block-SCC"},
            ChartLabels -> {labels, None},
            BarOrigin  -> Left,
            Frame -> True,
            FrameLabel -> {"speedup vs Lexicographic", ""},
            GridLines -> {{1.0}, None},
            ImageSize -> {500, 800},
            PlotLabel -> "Block-ordering speedup (>1.0 = faster than Lexicographic)"]
    ];

speedupBarChart


(* ::Section:: *)
(*Section 3 \[LongDash] Per-scenario detail*)


(* ::Text:: *)
(*Each subsection below assembles the available data for one scenario (top-N branching disjuncts, structural-role rows, sensitivity rows, benchmark times). Cells are independent \[LongDash] open and evaluate only the scenarios you care about.*)


scenarioDetail[scen_String] :=
    Module[{br = branchingData[scen], st = structuralData[scen],
            sn = sensitivityData[scen], bench = scenarioWallTimes[scen],
            sections = {}},
        AppendTo[sections, Style["Scenario: " <> scen, 15, Bold]];
        If[AssociationQ[bench] && bench =!= <||>,
            AppendTo[sections,
                Column[{Style["Wall time per ordering (s)", 12, Bold],
                        Dataset[bench]}, Spacings -> 0.3]]];
        If[AssociationQ[br],
            AppendTo[sections,
                Column[{Style["Branching totals", 12, Bold],
                        Dataset[Lookup[br, "totals", <||>]],
                        Style["Top-10 substituted disjuncts", 12, Bold],
                        Dataset[Take[Lookup[br, "perDisjunct", {}],
                                     UpTo[10]]]}, Spacings -> 0.3]]];
        If[AssociationQ[st],
            AppendTo[sections,
                Column[{Style["Structural rows (per abstract complementarity)", 12, Bold],
                        Style["outcomes: " <>
                              ToString[Lookup[st, "outcomeCounts", <||>]], 11],
                        Dataset[Take[Lookup[st, "rows", {}], UpTo[20]]]},
                       Spacings -> 0.3]]];
        If[AssociationQ[sn] && !TrueQ[Lookup[sn, "skipped", False]],
            AppendTo[sections,
                Column[{Style["Sensitivity rows (per disjunct)", 12, Bold],
                        Dataset[Take[Lookup[sn, "rows", {}], UpTo[20]]]},
                       Spacings -> 0.3]],
            If[AssociationQ[sn],
                AppendTo[sections,
                    Style["Sensitivity skipped: " <>
                          ToString[Lookup[sn, "reason", "?"]], 11, Italic]]]];
        Column[sections, Spacings -> 1]
    ];


(* Convenience accessors. Evaluate scenarioDetail["case_22"] (or any key
   from Keys[branchingData]) to see one scenario in full. *)
allScenarioKeys := Sort @ DeleteDuplicates @ Join[
    Keys[branchingData], Keys[structuralData], Keys[sensitivityData]];

allScenarioKeys


(* Example: render the top-branching scenario. *)
scenarioDetail["case_22"]


(* ::Section:: *)
(*Section 4 \[LongDash] Decision*)


(* ::Text:: *)
(*Three forks were posed in the original plan. The empirical answer for each, based on PRs 1-3 data:*)


(* ::Text:: *)
(*Fork 1 \[LongDash] Strong topology rule (P(outcome | feature bin) > 0.95 at meaningful n). NOT MET. The strongest cached-only signal is P(Active | junctionDegree<=2) = 0.80 (n = 218); P(Active | boundaryDistance = 0) = 0.50 dropping to 0.20 at distance 3. Useful as ordering heuristics, not as a hard preprocessing rule that pins variables before the solver runs.*)


(* ::Text:: *)
(*Fork 2 \[LongDash] Strong block-ordering winner. MET. 21 of 37 scenarios show >5% speedup under at least one Block-*; 0 show >10% regression. Best wins: Achdou_2023_junction 4.95\[Times], case_3 2.43\[Times], case_20 1.78\[Times], Grid0303 1.65\[Times], Grid0404 -5s on a 30s solve. Block-Edge and Block-SCC are roughly tied; Block-Vertex is competitive on tiny scenarios but loses on the bigger cracked cases.*)


(* ::Text:: *)
(*Fork 3 \[LongDash] No rule visible. Not the outcome.*)


(* ::Text:: *)
(*Bonus finding from the sensitivity sweep (323 disjuncts \[Times] 646 forced branches on the small subset): every modified solve produced the same flow-variable values as the baseline. Disjuncts shape the *search*, not the *answer*. This explains why block ordering helps without changing semantics: locality surfaces redundancy earlier.*)


(* ::Text:: *)
(*Recommended next step: ship "DisjunctOrdering" -> "Block-Edge" (or "Block-SCC") as the new default in a separate PR. Gate it on a BENCHMARKS.md refresh with both before/after numbers, since the orderings are opt-in today and a default flip is user-visible. The hard cases (Grid >= 5*5, HRF, case_21, case_23) are not cracked by ordering alone and remain on the oracle / numeric-fictitious-play track.*)
