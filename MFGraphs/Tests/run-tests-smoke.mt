(* Wolfram Language Test file *)
(* Smoke test for the curated test runner entrypoint. *)

Test[
    Module[{proc, script},
        script = "/Users/ribeirrd/Documents/GitHub/MFGraphs/Scripts/RunTests.wls";
        proc = RunProcess[{"wolframscript", "-file", script, "fast"}];
        proc["ExitCode"] === 0 && StringContainsQ[proc["StandardOutput"], "Suite: fast"]
    ]
    ,
    True
    ,
    TestID -> "RunTests.wls: fast suite smoke test"
]
