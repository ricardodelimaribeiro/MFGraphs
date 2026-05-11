#!/usr/bin/env wolframscript
(* CompareSingleCase.wls — run the five solver columns on ONE case in a fresh
   kernel. Designed for big grids (4x4, 4x5) where the multi-case driver
   accumulates enough session age to drop its MathLink connection.

   Usage:
     wolframscript -file Scripts/CompareSingleCase.wls grid-4x4
     wolframscript -file Scripts/CompareSingleCase.wls grid-4x5

   Recommended on macOS to suppress App Nap / idle sleep:
     caffeinate -dimsu wolframscript -file Scripts/CompareSingleCase.wls grid-4x4

   Output: prints one markdown table row to stdout, in the same column shape
   as Results/dnf_recursive_vs_procedural.md. Append it manually if desired. *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength = 0;

(* Per-call timeouts identical to CompareDNFRecursiveProcedural.wls. *)
$Timeout            = 300;
$BoolConvertTimeout = 60;
$BMTimeout          = 60;
$BMRTimeout         = 600;

$KnownCases = <|
    "grid-2x2" -> Function[gridScenario[{2, 2}, {{1, 100}}, {{4, 0}}]],
    "grid-3x2" -> Function[gridScenario[{3, 2}, {{1, 100}}, {{6, 0}}]],
    "grid-2x3" -> Function[gridScenario[{2, 3}, {{1, 100}}, {{6, 0}}]],
    "grid-2x4" -> Function[gridScenario[{2, 4}, {{1, 100}}, {{8, 0}}]],
    "grid-3x3" -> Function[gridScenario[{3, 3}, {{1, 100}}, {{9, 0}}]],
    "grid-3x4" -> Function[gridScenario[{3, 4}, {{1, 100}}, {{12, 0}}]],
    "grid-4x4" -> Function[gridScenario[{4, 4}, {{1, 100}}, {{16, 0}}]],
    "grid-4x5" -> Function[gridScenario[{4, 5}, {{1, 100}}, {{20, 0}}]]
|>;

args = Rest[$ScriptCommandLine];
If[Length[args] < 1,
    Print["Usage: wolframscript -file Scripts/CompareSingleCase.wls <case-name>"];
    Print["Known cases: ", Keys[$KnownCases]];
    Exit[1]
];
caseName = First[args];
If[!KeyExistsQ[$KnownCases, caseName],
    Print["Unknown case '", caseName, "'. Known: ", Keys[$KnownCases]];
    Exit[1]
];

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

fmtNum[x_?NumericQ] := ToString[NumberForm[N[x], {Infinity, 3}]];
fmtNum[x_]          := ToString[x];

builder = $KnownCases[caseName];
Print["[", caseName, "] building scen..."];
scen = builder[];
sys      = makeSystem[scen];
inputs   = solversTools`Private`buildSolverInputs[sys];
{constraints, allVars} = inputs[[1 ;; 2]];
nVars     = Length[allVars];
nConj     = If[Head[constraints] === And, Length[constraints], 1];
Print["  vars: ", nVars, "  conjuncts: ", nConj];

(* dnfReduce (recursive) *)
{tRec, resRec} = AbsoluteTiming[
    TimeConstrained[dnfReduce[True, constraints], $Timeout, $Aborted]
];
Print["  recursive:  ", tRec, "s",
      If[resRec === $Aborted, " (TIMEOUT)", ""]];

(* dnfReduceProcedural *)
{tProc, resProc} = AbsoluteTiming[
    TimeConstrained[
        solversTools`Private`dnfReduceProcedural[True, constraints, allVars],
        $Timeout, $Aborted
    ]
];
Print["  procedural: ", tProc, "s",
      If[resProc === $Aborted, " (TIMEOUT)", ""]];

recCan       = If[resRec  === $Aborted, $Aborted, canonicalBranches[resRec]];
procCan      = If[resProc === $Aborted, $Aborted, canonicalBranches[resProc]];
recBranches  = If[resRec  === $Aborted, "?", Length[recCan]];
procBranches = If[resProc === $Aborted, "?", Length[procCan]];
match = Which[
    resRec === $Aborted || resProc === $Aborted, "TIMEOUT",
    recCan === procCan,                          "OK",
    True,                                        "MISMATCH"
];
Print["  verdict: ", match];

(* Free large intermediate results before the Boolean stage. On grid-4x4 / 4x5
   the dnfReduce outputs (7200 / 27000 branches) plus their canonical sorted
   forms hold enough memory that the kernel can be OOM-killed during the
   subsequent BooleanConvert allocation. wolframscript reports an OOM kernel
   loss as "license error", which is misleading. *)
ClearAll[resRec, resProc, recCan, procCan];
Share[];
Print["  freed dnf results; mem now ", MemoryInUse[]];

(* BooleanConvert (Stage 1 only) *)
{tBool, resBool} = AbsoluteTiming[
    TimeConstrained[
        BooleanConvert[constraints, "DNF"],
        $BoolConvertTimeout, $Aborted
    ]
];
boolVerdict = If[resBool === $Aborted, "TIMEOUT", "OK"];
boolBranches = Which[
    resBool === $Aborted, "?",
    Head[resBool] === Or, ToString[Length[resBool]],
    True,                 "1"
];
Print["  BooleanConvert: ", tBool, "s",
      If[resBool === $Aborted, " (TIMEOUT)", ""],
      "  disjuncts: ", boolBranches];

(* booleanMinimizeSystem *)
{tBM, resBM} = AbsoluteTiming[
    TimeConstrained[
        Quiet[booleanMinimizeSystem[sys, "DisjunctTimeout" -> 10],
              {booleanMinimizeSystem::multisol}],
        $BMTimeout, $Aborted
    ]
];
bmVerdict = Which[
    resBM === $Aborted,                  "TIMEOUT",
    ListQ[resBM] || AssociationQ[resBM], "OK",
    True,                                "FAIL"
];
Print["  BooleanMinimize: ", tBM, "s  verdict: ", bmVerdict];

(* booleanMinimizeReduceSystem — runs in a fresh kernel so that OOM during BMR
   cannot kill the current kernel and lose the rows already computed. *)
ClearAll[resBool]; Share[];
bmrScript = FileNameJoin[{$RepoRoot, "Scripts", "CompareSingleCaseBMR.wls"}];
bmrProc = RunProcess[{"wolframscript", "-file", bmrScript, caseName}];
bmrOut  = Lookup[bmrProc, "StandardOutput", ""];
bmrLine = SelectFirst[StringSplit[bmrOut, "\n"], StringStartsQ[#, "BMR_RESULT|"] &, ""];
{bmrVerdict, tBMR} = If[bmrLine === "",
    {"FAIL", 0.0},
    Module[{parts = StringSplit[bmrLine, "|"]},
        {parts[[2]], ToExpression[parts[[3]]]}
    ]
];
Print["  BooleanMinimizeReduce: ", tBMR, "s  verdict: ", bmrVerdict, "  (subprocess)"];

(* Single markdown row, same column shape as the multi-case report. *)
row = StringJoin[
    "| ", caseName,
    " | ", ToString[nVars],
    " | ", match,
    " | ", ToString[recBranches],
    " | ", ToString[procBranches],
    " | ", fmtNum[tRec],
    " | ", fmtNum[tProc],
    " | ", boolVerdict,
    " | ", boolBranches,
    " | ", fmtNum[tBool],
    " | ", bmVerdict,
    " | ", fmtNum[tBM],
    " | ", bmrVerdict,
    " | ", fmtNum[tBMR],
    " |"
];
Print[];
Print["MARKDOWN ROW:"];
Print[row];

Exit[
    If[match === "MISMATCH", 1, 0]
];
