#!/usr/bin/env wolframscript
(* ProfileEquivCheck.wls -- Profile equivalence check time for a single case *)
(* Usage: wolframscript -file Scripts/ProfileEquivCheck.wls <case_key> [timeout_seconds] *)
(* Keys: Attraction_noswitch, Attraction_switch, Braess *)

$caseKey = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], "Attraction_noswitch"];
$equivTimeout = If[Length[$ScriptCommandLine] >= 3, ToExpression[$ScriptCommandLine[[3]]], 600];

Print["=== Equivalence Check Profiler ==="];
Print["Case: ", $caseKey];
Print["Timeout: ", $equivTimeout, "s"];
Print["Date: ", DateString[]];
Print[];

PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

Get[FileNameJoin[{Directory[], "Scripts", "BenchmarkHelpers.wls"}]];

$ComparisonCases = <|
  "Attraction_noswitch" -> {12, {I1 -> 100, U1 -> 0}},
  "Attraction_switch" -> {11, {I1 -> 100, U1 -> 0, S1 -> 1, S2 -> 1, S3 -> 1, S4 -> 1}},
  "Braess" -> {"Braess congest", {}}
|>;

(* --- Extract DNF input (same as CompareDNF.wls) --- *)
ExtractDNFInput[eqs_Association] :=
  Module[{preEqs, newSystem, initRules, sorted, xp},
    preEqs = If[KeyExistsQ[eqs, "InitRules"], eqs, MFGPreprocessing[eqs]];
    initRules = Lookup[preEqs, "InitRules", <||>];
    newSystem = Lookup[preEqs, "NewSystem", {True, True, True}];
    Module[{js, approxJs, costpluscurrents, Ncpc},
      js = Lookup[preEqs, "js", {}];
      approxJs = AssociationThread[js, 0 js];
      costpluscurrents = Lookup[preEqs, "costpluscurrents", {}];
      Ncpc = Expand /@ (costpluscurrents /. approxJs);
      initRules = Expand /@ (initRules /. Ncpc);
      newSystem = newSystem /. Ncpc;
    ];
    {newSystem, initRules} = TripleClean[{newSystem, initRules}];
    If[newSystem[[3]] === True,
      {True, True, initRules, "trivial"},
      sorted = DeduplicateByComplexity[newSystem[[3]]];
      xp = And @@ Most[newSystem];
      {xp, sorted, initRules, "nontrivial"}
    ]
  ];

(* --- Main --- *)
$caseSpec = $ComparisonCases[$caseKey];
If[!ListQ[$caseSpec] || Length[$caseSpec] =!= 2,
  Print["ERROR: Unknown case key '", $caseKey, "'. Available keys: ", Keys[$ComparisonCases]];
  Quit[1]
];

$exampleKey = $caseSpec[[1]];
$params = $caseSpec[[2]];

Print["Loading case ", $exampleKey, "..."];
$data = GetExampleData[$exampleKey];
If[$params =!= {}, $data = $data /. $params];
$d2e = DataToEquations[$data];
ClearSolveCache[];

{$xp, $sorted, $initRules, $tag} = Quiet @ Check[ExtractDNFInput[$d2e], {$Failed, $Failed, <||>, "failed"}];
If[$tag === "trivial" || $tag === "failed",
  Print["  Case is ", $tag, ". Nothing to profile."];
  Quit[]
];

Print["  xp: ", LeafCount[$xp], " leaves, sorted: ", If[Head[$sorted] === Or, Length[$sorted], 1], " disjuncts"];

(* Run BooleanConvert *)
Print["  Running BooleanConvert..."];
{$bcTime, $bcResult} = AbsoluteTiming[
  TimeConstrained[Quiet @ Check[BooleanConvert[$xp && $sorted, "DNF"], $Failed], 300, $TimedOut]
];
Print["  BooleanConvert: ", If[$bcResult === $TimedOut, "TIMED OUT",
  ToString@N[$bcTime, 4] <> "s, " <> ToString[LeafCount[$bcResult]] <> " leaves"]];

(* Run DNFReduce *)
Print["  Running DNFReduce..."];
ClearSolveCache[];
{$dnfTime, $dnfResult} = AbsoluteTiming[
  TimeConstrained[Quiet @ Check[DNFReduce[$xp, $sorted], $Failed], 300, $TimedOut]
];
Print["  DNFReduce: ", If[$dnfResult === $TimedOut, "TIMED OUT",
  ToString@N[$dnfTime, 4] <> "s, " <> ToString[LeafCount[$dnfResult]] <> " leaves"]];

(* --- The main event: equivalence check with extended timeout --- *)
If[$bcResult === $TimedOut || $bcResult === $Failed || $dnfResult === $TimedOut || $dnfResult === $Failed,
  Print["\n  SKIPPING equivalence check (a method timed out or failed)."];
  Quit[]
];

Print["\n  Running Reduce[bc \\[Equivalent] dnf, Reals] with ", $equivTimeout, "s timeout..."];
{$equivTime, $equiv} = AbsoluteTiming[
  TimeConstrained[
    Quiet @ Check[Reduce[$bcResult \[Equivalent] $dnfResult, Reals], "ERROR"],
    $equivTimeout, "EQUIV_TIMEOUT"
  ]
];

Print["  Result: ", Which[
  $equiv === True, "EQUIVALENT",
  $equiv === "EQUIV_TIMEOUT", "EQUIV_TIMEOUT (still running at " <> ToString[$equivTimeout] <> "s)",
  $equiv === "ERROR", "EQUIV_ERROR",
  True, "NOT_EQUIVALENT: " <> ToString[InputForm[$equiv]]
]];
Print["  Time: ", ToString@N[$equivTime, 6], "s"];
