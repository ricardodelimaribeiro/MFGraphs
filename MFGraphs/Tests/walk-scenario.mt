(* The gate-by-gate scenario walk as a fast-suite regression. The /scenario
   skill and Scripts/WalkScenario.wls both run Scripts/WalkScenarioHelpers.wls,
   which calls the pipeline through public symbols only. If one of those
   symbols is renamed or changes its contract, a test here fails, so the skill
   cannot go stale unnoticed the way its first version did. *)

Module[{root},
    root = Which[
        StringQ[$RepoRoot], $RepoRoot,
        StringQ[$InputFileName] && $InputFileName =!= "",
            FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}],
        True, Directory[]
    ];
    Get[FileNameJoin[{root, "Scripts", "WalkScenarioHelpers.wls"}]]
];

walkStatuses[walk_] := #["Status"] & /@ walk["Gates"];

(* A fully determined example passes every gate. *)
$walkCase12 = walkScenario[getExampleScenario[12, {{1, 80}}, {{4, 0}}]];

Test[
    {$walkCase12["Result"], walkStatuses[$walkCase12]}
    ,
    {"OK", {"PASS", "PASS", "PASS", "PASS", "PASS", "PASS"}}
    ,
    TestID -> "Walk: case 12 passes all six gates"
]

Test[
    #["Gate"] & /@ $walkCase12["Gates"]
    ,
    {"scenario", "unknowns", "system", "solve", "valid", "exact"}
    ,
    TestID -> "Walk: gates are reported in pipeline order"
]

(* Completeness is searched only on request, and the report says so. *)
Test[
    {$walkCase12["Exact"]["Soundness"], $walkCase12["Exact"]["Completeness"]}
    ,
    {"Proved", "NotChecked"}
    ,
    TestID -> "Walk: exact gate proves soundness and leaves completeness NotChecked"
]

(* An underdetermined example is OK and keeps its residual. *)
$walkGrid0303 = walkScenario[getExampleScenario["Grid0303", {{1, 100}}, {{9, 0}}]];

Test[
    {$walkGrid0303["Result"],
     AssociationQ[$walkGrid0303["Solution"]] &&
        !MemberQ[{True, False}, Lookup[$walkGrid0303["Solution"], "Residual", True]]}
    ,
    {"OK", True}
    ,
    TestID -> "Walk: Grid0303 is OK with a residual present"
]

(* A raw association goes through makeScenario; an entry vertex that is also
   an exit vertex is rejected at the first gate and nothing later runs. *)
$walkInvalid = walkScenario[<|"Model" -> <|
    "Vertices" -> {1, 2}, "Adjacency" -> {{0, 1}, {0, 0}},
    "Entries" -> {{1, 10}}, "Exits" -> {{1, 0}}, "Switching" -> {}|>|>];

Test[
    {$walkInvalid["Result"], walkStatuses[$walkInvalid]}
    ,
    {"GateFailed", {"FAIL", "SKIP", "SKIP", "SKIP", "SKIP", "SKIP"}}
    ,
    TestID -> "Walk: an invalid model fails the scenario gate and skips the rest"
]

(* Non-critical congestion stops at the system gate, before any solver runs. *)
$walkAlpha2 = walkScenario[getExampleScenario[12, {{1, 80}}, {{4, 0}}, {}, 2]];

Test[
    {$walkAlpha2["Result"], walkStatuses[$walkAlpha2]}
    ,
    {"GateFailed", {"PASS", "PASS", "FAIL", "SKIP", "SKIP", "SKIP"}}
    ,
    TestID -> "Walk: Alpha = 2 fails the system gate"
]

(* A chain whose only through transition is blocked has no solution. That is
   a result of its own, not a failed gate. *)
$walkBlocked = walkScenario[gridScenario[{3}, {{1, 10}}, {{3, 0}}, {{1, 2, 3, Infinity}}]];

Test[
    {$walkBlocked["Result"], walkStatuses[$walkBlocked]}
    ,
    {"NoSolution", {"PASS", "PASS", "PASS", "NOSOLUTION", "SKIP", "SKIP"}}
    ,
    TestID -> "Walk: a blocked chain reports NoSolution at the solve gate"
]

(* The text report has one line per gate and ends on the result. *)
Test[
    With[{lines = StringSplit[formatWalk[$walkCase12], "\n"]},
        {Length[lines], StringStartsQ[First[lines], "[1] scenario"], Last[lines]}]
    ,
    {7, True, "RESULT        OK"}
    ,
    TestID -> "Walk: formatWalk prints six gate lines and a RESULT line"
]
