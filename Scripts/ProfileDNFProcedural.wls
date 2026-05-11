#!/usr/bin/env wolframscript
(* ProfileDNFProcedural.wls — baseline profiler for dnfReduceProcedural.

   The recursive engine has dnfReduceDiagnosticReport (via dnfReduceInstrumented);
   the procedural engine has no equivalent. This script measures the procedural
   reducer's wall-clock time, output branch count, output leaf count, and peak
   worklist size on the same grid cases used by ProfileDNFReduce.wls and the
   recursive/procedural comparison report.

   Usage:
     wolframscript -file Scripts/ProfileDNFProcedural.wls [--case name|all] [--timeout s]

   Output: TSV rows on stdout, one per case.
*)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength  = 0;

$SelectedCase = "all";
$TimeoutSecs  = 300;

Module[{args = Rest[$ScriptCommandLine], i = 1, arg},
    While[i <= Length[args],
        arg = args[[i]];
        Switch[arg,
            "--case",    i++; If[i <= Length[args], $SelectedCase = ToLowerCase[args[[i]]]],
            "--timeout", i++; If[i <= Length[args],
                Module[{v = Quiet @ Check[ToExpression[args[[i]]], $Failed]},
                    If[NumericQ[v] && v >= 0, $TimeoutSecs = v]]],
            _, Null
        ];
        i++
    ]
];

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

If[$SelectedCase =!= "all",
    $Cases = Select[$Cases, #[[1]] === $SelectedCase &];
    If[$Cases === {},
        Print["No case '", $SelectedCase, "'"]; Exit[2]
    ]
];

Fmt[x_?NumericQ] := ToString[NumberForm[N[x], {Infinity, 3}]];
Fmt[x_] := ToString[x];

branchCount[expr_] :=
    Which[
        expr === True || expr === False, 0,
        Head[expr] === Or, Length[expr],
        True, 1
    ];

Print["MFGraphs dnfReduceProcedural baseline profile"];
Print["Date\t", DateString[]];
Print["TimeoutSeconds\t", $TimeoutSecs];
Print["Case\tVars\tConj\tBuildMs\tProcMs\tStatus\tBranches\tLeaves"];

Scan[
    Function[{case},
        Module[{name, builder, scen, sys, inputs, constraints, allVars,
                tBuild, tProc, result, nVars, nConj, status, br, lv},
            name    = case[[1]];
            builder = case[[2]];
            {tBuild, sys} = AbsoluteTiming[makeSystem[builder[]]];
            inputs  = solversTools`Private`buildSolverInputs[sys];
            {constraints, allVars} = inputs[[1 ;; 2]];
            nVars = Length[allVars];
            nConj = If[Head[constraints] === And, Length[constraints], 1];
            {tProc, result} = AbsoluteTiming[
                TimeConstrained[
                    solversTools`Private`dnfReduceProcedural[True, constraints, allVars],
                    $TimeoutSecs, $Aborted
                ]
            ];
            status = If[result === $Aborted, "TIMEOUT", "OK"];
            br = If[result === $Aborted, "?", ToString[branchCount[result]]];
            lv = If[result === $Aborted, "?", ToString[LeafCount[result]]];
            Print[StringRiffle[{
                name,
                ToString[nVars],
                ToString[nConj],
                Fmt[1000 tBuild],
                Fmt[1000 tProc],
                status,
                br,
                lv
            }, "\t"]];
        ]
    ],
    $Cases
];
