#!/usr/bin/env wolframscript

PrependTo[$Path, Directory[]];

testsRoot = FileNameJoin[{Directory[], "MFGraphs", "Tests"}];
runTestsPath = FileNameJoin[{Directory[], "Scripts", "RunTests.wls"}];
runTestsContent = Import[runTestsPath, "Text"];

forbidden = {
    "NonLinearSolver",
    "MonotoneSolver",
    "MonotoneSolverFromData",
    "ClassifyKKT",
    "NormalizeEdgeFunction"
};

fastMatches = StringCases[
    runTestsContent,
    RegularExpression["\"([^\"]+\\.mt)\""] :> "$1"
];
fastFiles = DeleteDuplicates @ Select[fastMatches, StringEndsQ[#, ".mt"] &];

archiveLeak = Select[fastFiles, StringStartsQ[#, "archive/"] || StringContainsQ[#, "/archive/"] &];

activeTestPaths = Select[
    FileNames["*.mt", testsRoot, Infinity],
    !StringContainsQ[#, "/archive/"] &
];

symbolLeaks = Reap[
    Do[
        content = Import[path, "Text"];
        Do[
            If[StringContainsQ[content, sym],
                Sow[path <> " -> " <> sym]
            ],
            {sym, forbidden}
        ],
        {path, activeTestPaths}
    ]
][[2]];

symbolLeaks = If[symbolLeaks === {}, {}, First[symbolLeaks]];

If[archiveLeak === {} && symbolLeaks === {},
    Print["Critical-surface test guard: PASS"];
    Exit[0]
];

If[archiveLeak =!= {},
    Print["Critical-surface test guard: fast suite includes archived tests:"];
    Scan[Print["  - ", #] &, archiveLeak];
];

If[symbolLeaks =!= {},
    Print["Critical-surface test guard: active tests reference unsupported symbols:"];
    Scan[Print["  - ", #] &, symbolLeaks];
];

Exit[1];
