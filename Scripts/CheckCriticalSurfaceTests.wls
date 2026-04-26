#!/usr/bin/env wolframscript

$RepoRoot = "/Users/ribeirrd/Documents/GitHub/MFGraphs";

PrependTo[$Path, $RepoRoot];

testsRoot = FileNameJoin[{$RepoRoot, "MFGraphs", "Tests"}];
runTestsPath = FileNameJoin[{$RepoRoot, "Scripts", "RunTests.wls"}];
runTestsContent = Import[runTestsPath, "Text"];

forbidden = {
    "NonLinearSolver",
    "MonotoneSolver",
    "MonotoneSolverFromData",
    "ClassifyKKT",
    "NormalizeEdgeFunction"
};

ExtractSuiteFiles[content_String, suiteName_String] := Module[
    {suiteBlock, matches},
    suiteBlock = FirstCase[
        StringCases[
            content,
            RegularExpression["(?s)\"" <> suiteName <> "\"\\s*->\\s*\\{(.*?)\\}"] :> "$1"
        ],
        _String,
        Missing["SuiteNotFound", suiteName]
    ];
    If[MatchQ[suiteBlock, _Missing], Return[suiteBlock]];
    matches = StringCases[suiteBlock, RegularExpression["\"([^\"]+\\.mt)\""] :> "$1"];
    DeleteDuplicates @ Select[matches, StringEndsQ[#, ".mt"] &]
];

fastFiles = ExtractSuiteFiles[runTestsContent, "fast"];
If[MatchQ[fastFiles, _Missing],
    Print["Critical-surface test guard: could not locate fast suite declaration in Scripts/RunTests.wls"];
    Exit[1];
];

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
