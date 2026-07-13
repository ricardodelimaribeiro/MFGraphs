(* Usage-arity lint as a fast-suite regression (issue #242): a definition
   arity change without a docs regeneration must fail here, not lie dormant
   until the next manual GenerateDocs.wls run. Shares its implementation
   with Scripts/GenerateDocs.wls via Scripts/UsageArityLint.wls. *)

Module[{root},
    root = Which[
        StringQ[$RepoRoot], $RepoRoot,
        StringQ[$InputFileName] && $InputFileName =!= "",
            FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}],
        True, Directory[]
    ];
    Get[FileNameJoin[{root, "Scripts", "UsageArityLint.wls"}]]
];

$usageArityReport = usageArityReport[];

(* Expected value {} means a failure prints the offending mismatch
   descriptions directly in the test report. *)
Test[
    $usageArityReport["Issues"]
    ,
    {}
    ,
    TestID -> "Usage arity: every checkable documented call pattern matches its definition"
]

Test[
    $usageArityReport["Warnings"]
    ,
    {}
    ,
    TestID -> "Usage arity: no unterminated call groups in usage strings"
]

(* Decay guard: if the definition introspection stops matching the package's
   definition style, the lint silently checks nothing -- fail loudly instead.
   82 symbols are checked as of 2026-07; the floor is deliberately loose. *)
Test[
    $usageArityReport["Checked"] >= 50
    ,
    True
    ,
    TestID -> "Usage arity: coverage has not decayed (at least 50 symbols checked)"
]
