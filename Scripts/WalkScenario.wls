#!/usr/bin/env wolframscript
(* Usage:
     wolframscript -file Scripts/WalkScenario.wls --list
     wolframscript -file Scripts/WalkScenario.wls --case <key> [--solve-timeout s] [--completeness s]
     wolframscript -file Scripts/WalkScenario.wls --scenario-file <path.wl> [--solve-timeout s] [--completeness s]

   Walks one scenario through the solver pipeline gate by gate (scenario,
   unknowns, system, solve, valid, exact), prints one line per gate, and stops
   at the first gate that does not pass. The walk itself lives in
   Scripts/WalkScenarioHelpers.wls. This is the runner behind the /scenario skill.

   --list                    print the built-in example keys
   --case <key>              a built-in example with its default entries and exits
                             ($ExampleBoundaryDefaults in Scripts/SolutionCacheHelpers.wls)
   --scenario-file <path.wl> a file whose value is a scenario or a raw association
                             for makeScenario; it is read after MFGraphs` loads, so
                             gridScenario, getExampleScenario etc. work inside it
   --solve-timeout s         seconds allowed for the solve gate (default 60)
   --completeness s          also run the exact completeness search, s seconds at
                             most (off by default; its status is printed verbatim)

   Exit codes: 0 all gates passed, 1 a gate failed, 2 usage error,
               3 no solution, 4 solve timeout.

   The flags are deliberately not named --timeout or --file. wolframscript
   acts on a --timeout argument itself, even when it follows the script name:
   it stops the whole run after that many seconds, prints a stray token and
   exits 0, so the solve gate could never report its own timeout. *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
Get[FileNameJoin[{$RepoRoot, "Scripts", "SolutionCacheHelpers.wls"}]];
Get[FileNameJoin[{$RepoRoot, "Scripts", "WalkScenarioHelpers.wls"}]];

$UsageLine = "Usage: wolframscript -file Scripts/WalkScenario.wls (--list | --case <key> | --scenario-file <path.wl>) [--solve-timeout s] [--completeness s]";

UsageError[message_String] := (Print[message]; Print[$UsageLine]; Exit[2]);

PositiveNumberArg[flag_String, text_String] :=
  With[{value = Quiet @ Check[ToExpression[text], $Failed]},
    If[NumericQ[value] && value > 0, value,
      UsageError[flag <> " needs a positive number, got: " <> text]]
  ];

$ListOnly = False;
$CaseKey = None;
$ScenarioFile = None;
$SolveTimeoutSecs = 60;
$CompletenessSecs = None;

Module[{args = Rest[$ScriptCommandLine], i = 1, arg, next},
  next[flag_] := (i++; If[i <= Length[args], args[[i]], UsageError[flag <> " needs a value"]]);
  While[i <= Length[args],
    arg = args[[i]];
    Which[
      arg === "--list", $ListOnly = True,
      arg === "--case", $CaseKey = next["--case"],
      arg === "--scenario-file", $ScenarioFile = next["--scenario-file"],
      arg === "--solve-timeout", $SolveTimeoutSecs = PositiveNumberArg["--solve-timeout", next["--solve-timeout"]],
      arg === "--completeness", $CompletenessSecs = PositiveNumberArg["--completeness", next["--completeness"]],
      StringEndsQ[arg, ".wls"], Null,
      True, UsageError["Unknown argument: " <> arg]
    ];
    i++
  ]
];

If[$ListOnly,
  Scan[Print, listExampleScenarios[]];
  Exit[0]
];

If[($CaseKey === None) === ($ScenarioFile === None),
  UsageError["Give exactly one of --case or --scenario-file."]
];

If[$CaseKey =!= None,
  (* Registry keys are integers or strings; "12" on the command line means 12. *)
  $CaseKey = If[StringMatchQ[$CaseKey, DigitCharacter ..], ToExpression[$CaseKey], $CaseKey];
  If[!MemberQ[listExampleScenarios[], $CaseKey],
    UsageError["No built-in example named '" <> ToString[$CaseKey] <> "'. Run with --list to see the keys."]
  ];
  $WalkLabel = ToString[$CaseKey];
  $WalkInput = getExampleScenario[$CaseKey, Sequence @@ $ExampleBoundaryDefaults[$CaseKey]]
  ,
  If[!FileExistsQ[$ScenarioFile],
    UsageError["Scenario file not found: " <> $ScenarioFile]
  ];
  $WalkLabel = $ScenarioFile;
  $WalkInput = Get[$ScenarioFile]
];

(* The header goes out before the walk starts. If the kernel terminates the
   evaluation mid-walk (it can still exceed the raised recursion limit, or
   crash), the output ends without a RESULT line while the exit code stays 0,
   so a missing RESULT line is the signal, not the exit code. *)
Print["SCENARIO      ", $WalkLabel];

$Walk = walkScenario[$WalkInput, "Timeout" -> $SolveTimeoutSecs, "Completeness" -> $CompletenessSecs];

Print[formatWalk[$Walk]];

Exit[
  Switch[$Walk["Result"],
    "OK", 0,
    "NoSolution", 3,
    "Timeout", 4,
    _, 1
  ]
];

(* :!CodeAnalysis::EndBlock:: *)
