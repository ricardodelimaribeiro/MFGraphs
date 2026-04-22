#!/usr/bin/env wolframscript
(* LEGACY_GETEXAMPLEDATA: This script still depends on GetExampleData and is excluded from the active scenario-first workflow. *)
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

$RepoRoot = "/Users/ribeirrd/Documents/GitHub/MFGraphs";
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;

Get[FileNameJoin[{$RepoRoot, "Scripts", "BenchmarkHelpers.wls"}]];

(* --- Main --- *)
$caseSpec = Lookup[$NamedBenchmarkCases, $caseKey, $Failed];
If[$caseSpec === $Failed,
  Print["ERROR: Unknown case key '", $caseKey, "'. Available keys: ", Keys[$NamedBenchmarkCases]];
  Quit[1]
];

$exampleKey = $caseSpec[[1]];
$params = $caseSpec[[2]];

Print["Loading case ", $exampleKey, "..."];
$data = GetExampleData[$exampleKey];
If[$params =!= {}, $data = $data /. $params];
$d2e = DataToEquations[$data];
ClearSolveCache[];

(* Extract DNF input using centralized helper *)
{$xp, $sorted, $initRules, $tag} = Quiet @ Check[ExtractDNFInput[$d2e], {$Failed, $Failed, <||>, "failed"}];
If[$tag === "trivial" || $tag === "failed",
  Print["  Case is ", $tag, ". Nothing to profile."];
  Quit[]
];

Print["  xp: ", LeafCount[$xp], " leaves, sorted: ", CountDisjuncts[$sorted], " disjuncts"];

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
    $equivTimeout,
    "TIMEOUT"
  ]
];

Print["  Equivalence: ", $equiv, " (", ToString@N[$equivTime, 4], "s)"];
If[$equiv === "TIMEOUT",
  Print["  [FAIL] Equivalence check exceeded the extended timeout of ", $equivTimeout, "s."]
];
If[$equiv === False,
  Print["  [FAIL] DNFReduce is NOT EQUIVALENT to BooleanConvert for this case!"]
];

Print["\nProfile completed at ", DateString[]];
