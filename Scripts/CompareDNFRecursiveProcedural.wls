#!/usr/bin/env wolframscript
(* Compare correctness of dnfReduce (recursive) vs dnfReduceProcedural (stack-based)
   on the real MFG cases used by ProfileDNFReduce.wls. Writes a markdown report to
   Results/dnf_recursive_vs_procedural.md with a per-case verdict.

   Equivalence rule: convert each output to a canonical sorted set of And-branches
   (order-insensitive within each branch and across branches). If canonical forms
   match, the two reducers produced the same DNF up to ordering. *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength = 0;

$Timeout = 300;

$Cases = {
    {"grid-2x2", Function[gridScenario[{2, 2}, {{1, 100}}, {{4, 0}}]]},
    {"grid-3x2", Function[gridScenario[{3, 2}, {{1, 100}}, {{6, 0}}]]},
    {"grid-2x3", Function[gridScenario[{2, 3}, {{1, 100}}, {{6, 0}}]]},
    {"grid-2x4", Function[gridScenario[{2, 4}, {{1, 100}}, {{8, 0}}]]},
    {"grid-3x3", Function[gridScenario[{3, 3}, {{1, 100}}, {{9, 0}}]]},
    {"grid-3x4", Function[gridScenario[{3, 4}, {{1, 100}}, {{12, 0}}]]},
    {"grid-4x4", Function[gridScenario[{4, 4}, {{1, 100}}, {{16, 0}}]]},
    {"grid-4x5", Function[gridScenario[{4, 5}, {{1, 100}}, {{20, 0}}]]}
};

(* --- canonicalize an Or-of-And DNF expression --- *)
canonicalBranches[True]  := {True};
canonicalBranches[False] := {};
canonicalBranches[expr_] := Module[{branches},
    branches = If[Head[expr] === Or, List @@ expr, {expr}];
    Sort[
        Map[
            Function[branch,
                Module[{conj},
                    conj = If[Head[branch] === And, List @@ branch, {branch}];
                    Sort[conj]
                ]
            ],
            branches
        ]
    ]
];

(* --- run a single case --- *)
runCase[name_String, builder_] := Module[
    {scenario, sys, inputs, constraints, allVars, tRec, resRec, tProc, resProc,
     match, recCan, procCan, recBranches, procBranches},
    Print["[", name, "] building scenario..."];
    scenario = builder[];
    sys = makeSystem[scenario];
    inputs = solversTools`Private`buildSolverInputs[sys];
    {constraints, allVars} = inputs[[1 ;; 2]];

    Print["  vars: ", Length[allVars], "  conjuncts: ",
          If[Head[constraints] === And, Length[constraints], 1]];

    {tRec, resRec} = AbsoluteTiming[
        TimeConstrained[
            dnfReduce[True, constraints],
            $Timeout, $Aborted
        ]
    ];
    Print["  recursive:  ", tRec, "s",
          If[resRec === $Aborted, " (TIMEOUT)", ""]];

    {tProc, resProc} = AbsoluteTiming[
        TimeConstrained[
            solversTools`Private`dnfReduceProcedural[True, constraints, allVars],
            $Timeout, $Aborted
        ]
    ];
    Print["  procedural: ", tProc, "s",
          If[resProc === $Aborted, " (TIMEOUT)", ""]];

    recCan       = If[resRec === $Aborted, $Aborted, canonicalBranches[resRec]];
    procCan      = If[resProc === $Aborted, $Aborted, canonicalBranches[resProc]];
    recBranches  = If[resRec === $Aborted, "?", Length[recCan]];
    procBranches = If[resProc === $Aborted, "?", Length[procCan]];
    match = Which[
        resRec === $Aborted || resProc === $Aborted, "TIMEOUT",
        recCan === procCan,                          "OK",
        True,                                        "MISMATCH"
    ];

    Print["  verdict: ", match];

    <|
        "Case"            -> name,
        "Vars"            -> Length[allVars],
        "Match"           -> match,
        "RecBranches"     -> recBranches,
        "ProcBranches"    -> procBranches,
        "RecSeconds"      -> N[tRec],
        "ProcSeconds"     -> N[tProc]
    |>
];

results = runCase @@@ $Cases;

(* --- write the markdown report --- *)
fmtNum[x_?NumericQ] := ToString[NumberForm[N[x], {Infinity, 3}]];
fmtNum[x_]          := ToString[x];

reportLines = Join[
    {
        "# DNF reducer correctness: recursive vs procedural",
        "",
        "Comparison of `dnfReduce` (recursive, public in `solversTools``) and `dnfReduceProcedural` (stack-based, private in `solversTools`Private``) on grid scenarios.",
        "",
        "Equivalence rule: each output is canonicalized to a sorted set of `And`-branches (order-insensitive within each branch and across branches). `OK` means both procedures produce the same DNF up to ordering.",
        "",
        "Timeout per call: " <> ToString[$Timeout] <> "s.",
        "",
        "| Case | Vars | Verdict | Recursive branches | Procedural branches | Recursive (s) | Procedural (s) |",
        "|---|---|---|---|---|---|---|"
    },
    Map[
        Function[r,
            StringJoin[
                "| ", r["Case"],
                " | ", ToString[r["Vars"]],
                " | ", r["Match"],
                " | ", ToString[r["RecBranches"]],
                " | ", ToString[r["ProcBranches"]],
                " | ", fmtNum[r["RecSeconds"]],
                " | ", fmtNum[r["ProcSeconds"]],
                " |"
            ]
        ],
        results
    ],
    {
        "",
        "Generated by `Scripts/CompareDNFRecursiveProcedural.wls`."
    }
];

outPath = FileNameJoin[{$RepoRoot, "Results", "dnf_recursive_vs_procedural.md"}];
If[!DirectoryQ[DirectoryName[outPath]], CreateDirectory[DirectoryName[outPath]]];
Export[outPath, StringRiffle[reportLines, "\n"], "Text"];
Print["Wrote ", outPath];

(* exit nonzero if any mismatch — useful for CI *)
If[MemberQ[results[[All, "Match"]], "MISMATCH"],
    Print["FAIL: at least one mismatch."]; Exit[1],
    Print["All cases match (or timed out)."]; Exit[0]
];
