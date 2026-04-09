#!/usr/bin/env wolframscript
(* JamaratFeasibilityPhaseDiagram.wls
   Export feasibility phase-diagram data for the Jamaratv9 case.

   Usage:
     wolframscript -file Scripts/JamaratFeasibilityPhaseDiagram.wls
     wolframscript -file Scripts/JamaratFeasibilityPhaseDiagram.wls profile=paper
     wolframscript -file Scripts/JamaratFeasibilityPhaseDiagram.wls profile=dense timeout=180

   Optional overrides (key=value):
     profile=quick|paper|dense
     timeout=<seconds>
     total=<I1+I2 total for flow-split sweep>
     basei1=<fixed I1 for cost sweeps>
     basei2=<fixed I2 for cost sweeps>
     baseu1=<baseline U1 for flow split + cost windows>
     baseu2=<baseline U2 for flow split + cost windows>
     baseu3=<baseline U3 for flow split + cost windows>
     anchoru1=<U1 anchor for cost-difference grid>
     flow=min:max:step
     d21=min:max:step
     d31=min:max:step
     u1=min:max:step
     u2=min:max:step
     u3=min:max:step

   Exports (Results/):
     jamarat_phase_points_<timestamp>.{csv,json}
     jamarat_phase_points_latest.{csv,json}
     jamarat_phase_boundaries_<timestamp>.{csv,json}
     jamarat_phase_boundaries_latest.{csv,json}
     jamarat_phase_costdiff_matrix_<timestamp>.csv
     jamarat_phase_costdiff_matrix_latest.csv
     jamarat_phase_metadata_<timestamp>.json
     jamarat_phase_metadata_latest.json
*)

Print["MFGraphs Jamarat Feasibility Phase Diagram Export"];
Print["Date: ", DateString[]];
Print["Mathematica: ", $Version];
Print[];

PrependTo[$Path, Directory[]];
Needs["MFGraphs`"];
Get[FileNameJoin[{Directory[], "Scripts", "BenchmarkHelpers.wls"}]];
MFGraphs`$MFGraphsVerbose = False;

ParseNumber[value_String] :=
    Module[{held},
        held =
            Quiet @ Check[
                ToExpression[StringTrim[value], InputForm, HoldComplete],
                $Failed
            ];
        Which[
            held === $Failed, $Failed,
            MatchQ[held, HoldComplete[_Integer | _Real | _Rational]],
                N @ ReleaseHold[held],
            True, $Failed
        ]
    ];
ParseNumber[value_?NumericQ] := N[value];
ParseNumber[_] := $Failed;

ParseRangeSpec[value_String] :=
    Module[{parts, parsed},
        parts = StringSplit[value, ":"];
        If[Length[parts] =!= 3, Return[$Failed, Module]];
        parsed = ParseNumber /@ parts;
        If[VectorQ[parsed, NumericQ] && parsed[[3]] =!= 0, parsed, $Failed]
    ];

ParseArgs[args_List] :=
    Module[{parsed = <||>, i = 1, arg, key, value, kv},
        While[i <= Length[args],
            arg = args[[i]];
            If[StringStartsQ[arg, "--"],
                key = ToLowerCase[StringDrop[arg, 2]];
                i++;
                If[i <= Length[args], parsed[key] = args[[i]]],
                If[StringContainsQ[arg, "="],
                    kv = StringSplit[arg, "=", 2];
                    If[Length[kv] == 2,
                        key = ToLowerCase[kv[[1]]];
                        value = kv[[2]];
                        parsed[key] = value
                    ]
                ]
            ];
            i++;
        ];
        parsed
    ];

RangeSpecToList[{min_, max_, step_}] := Range[min, max, step];

ApplyNumericOverride[config_Association, raw_Association, rawKey_String, configKey_String] :=
    Module[{parsed, updated},
        updated = config;
        If[KeyExistsQ[raw, rawKey],
            parsed = ParseNumber[raw[rawKey]];
            If[NumericQ[parsed],
                updated = Join[updated, <|configKey -> parsed|>],
                Print["Warning: invalid numeric value for ", rawKey, "=", raw[rawKey]]
            ]
        ];
        updated
    ];

ApplyRangeOverride[config_Association, raw_Association, rawKey_String, configKey_String] :=
    Module[{parsed, updated},
        updated = config;
        If[KeyExistsQ[raw, rawKey],
            parsed = ParseRangeSpec[raw[rawKey]];
            If[parsed =!= $Failed,
                updated = Join[updated, <|configKey -> parsed|>],
                Print["Warning: invalid range for ", rawKey, "=", raw[rawKey], " (expected min:max:step)"]
            ]
        ];
        updated
    ];

RunSolverWithTimeout[data_Association, timeout_?NumericQ] :=
    Module[{result, time, memBefore, memAfter, memDelta, status},
        memBefore = MemoryInUse[];
        {time, result} = AbsoluteTiming[
            TimeConstrained[
                Quiet[
                    MFGraphs`CriticalCongestionSolver[
                        MFGraphs`DataToEquations[data]
                    ],
                    {MFGraphs`MFGSystemSolver::nosolution}
                ],
                timeout,
                $TimedOut
            ]
        ];
        (* Per-run working-set delta (MaxMemoryUsed is session-global/monotonic). *)
        memAfter = MemoryInUse[];
        memDelta = Max[0, memAfter - memBefore];
        status = Which[
            result === $TimedOut, "TIMEOUT",
            AssociationQ[result], "OK",
            True, "FAILED"
        ];
        <|
            "Time" -> time,
            "Memory" -> memDelta,
            "Status" -> status,
            "Result" -> result
        |>
    ];

EvaluateScenario[scenario_Association, timeout_?NumericQ] :=
    Module[{data, wrapped, result, executionStatus, feasible, status, resultKind, message,
            i1, i2, u1, u2, u3},
        i1 = Lookup[scenario, "I1", 0];
        i2 = Lookup[scenario, "I2", 0];
        u1 = Lookup[scenario, "U1", 0];
        u2 = Lookup[scenario, "U2", 0];
        u3 = Lookup[scenario, "U3", 0];

        MFGraphs`ClearSolveCache[];
        data =
            MFGraphs`GetExampleData["Jamaratv9"] /. {
                MFGraphs`I1 -> i1,
                MFGraphs`I2 -> i2,
                MFGraphs`U1 -> u1,
                MFGraphs`U2 -> u2,
                MFGraphs`U3 -> u3
            };

        wrapped = RunSolverWithTimeout[data, timeout];

        result = Lookup[wrapped, "Result", Null];
        executionStatus = Lookup[wrapped, "Status", "UNKNOWN"];
        feasible = AssociationQ[result] && TrueQ[MFGraphs`IsFeasible[result]];
        status =
            Which[
                AssociationQ[result], Lookup[result, "Status", Lookup[result, "Feasibility", "Unknown"]],
                executionStatus === "TIMEOUT", "Timeout",
                executionStatus === "FAILED", "Failed",
                True, "Unknown"
            ];
        resultKind =
            If[AssociationQ[result],
                Lookup[result, "ResultKind", Missing["NotAvailable"]],
                Missing["NotAvailable"]
            ];
        message =
            If[AssociationQ[result],
                Lookup[result, "Message", None],
                Which[
                    executionStatus === "TIMEOUT", "Timeout",
                    executionStatus === "FAILED", "Failed",
                    True, None
                ]
            ];

        Join[
            scenario,
            <|
                "TotalInflow" -> i1 + i2,
                "FlowShareI1" -> If[i1 + i2 == 0, Indeterminate, N[i1/(i1 + i2)]],
                "Delta21" -> u2 - u1,
                "Delta31" -> u3 - u1,
                "Feasible" -> feasible,
                "Status" -> status,
                "ResultKind" -> resultKind,
                "Message" -> message,
                "ExecutionStatus" -> executionStatus,
                "SolveTime" -> Lookup[wrapped, "Time", Missing["NotAvailable"]],
                "SolveMemory" -> Lookup[wrapped, "Memory", Missing["NotAvailable"]]
            |>
        ]
    ];

EvaluateSlice[sliceLabel_String, scenarios_List, timeout_?NumericQ] :=
    Module[{n, checkpoint},
        n = Length[scenarios];
        Print["Running ", sliceLabel, " sweep (", n, " points)..."];
        checkpoint = Max[1, Ceiling[n/10]];
        MapIndexed[
            Function[{scenario, idx},
                If[idx[[1]] == 1 || idx[[1]] == n || Mod[idx[[1]], checkpoint] == 0,
                    Print["  ", sliceLabel, ": ", idx[[1]], "/", n]
                ];
                EvaluateScenario[scenario, timeout]
            ],
            scenarios
        ]
    ];

FeasibleIntervals[rows_List, variableKey_String] :=
    Module[{sorted, values, feasible, feasibleBlocks},
        If[rows === {}, Return[{}, Module]];
        sorted = SortBy[rows, Lookup[#, variableKey, 0] &];
        values = Lookup[sorted, variableKey, Missing["NotAvailable"]];
        feasible = Lookup[sorted, "Feasible", False];
        feasibleBlocks =
            Select[
                SplitBy[Range[Length[feasible]], feasible[[#]] &],
                TrueQ[feasible[[First[#]]]] &
            ];
        ({values[[First[#]]], values[[Last[#]]]} &) /@ feasibleBlocks
    ];

BoundaryRow[rows_List, slice_String, variableKey_String, fixed_Association] :=
    Module[{sorted, values, feasibleFlags, feasibleValues, infeasibleValues, intervals},
        sorted = SortBy[rows, Lookup[#, variableKey, 0] &];
        values = Lookup[sorted, variableKey, Missing["NotAvailable"]];
        feasibleFlags = Lookup[sorted, "Feasible", False];
        feasibleValues = Pick[values, feasibleFlags, True];
        infeasibleValues = Pick[values, feasibleFlags, False];
        intervals = FeasibleIntervals[rows, variableKey];
        Join[
            <|
                "Slice" -> slice,
                "Variable" -> variableKey,
                "SampledPoints" -> Length[rows],
                "FeasiblePoints" -> Count[feasibleFlags, True],
                "InfeasiblePoints" -> Count[feasibleFlags, False],
                "FirstFeasible" -> If[feasibleValues === {}, Missing["NotAvailable"], First[feasibleValues]],
                "LastFeasible" -> If[feasibleValues === {}, Missing["NotAvailable"], Last[feasibleValues]],
                "FirstInfeasible" -> If[infeasibleValues === {}, Missing["NotAvailable"], First[infeasibleValues]],
                "LastInfeasible" -> If[infeasibleValues === {}, Missing["NotAvailable"], Last[infeasibleValues]],
                "FeasibleIntervals" -> intervals,
                "FeasibleIntervalsText" -> ToString[InputForm[intervals]]
            |>,
            fixed
        ]
    ];

CostGridBoundaryRows[costRows_List] :=
    Module[{d31Values, rowsAtD31, intervals},
        d31Values = Sort @ DeleteDuplicates[Lookup[costRows, "Delta31", Missing["NotAvailable"]]];
        Table[
            rowsAtD31 = Select[costRows, Lookup[#, "Delta31", Missing["NotAvailable"]] === d31 &];
            intervals = FeasibleIntervals[rowsAtD31, "Delta21"];
            <|
                "Slice" -> "CostDifferenceGrid",
                "Variable" -> "Delta21",
                "Delta31" -> d31,
                "SampledPoints" -> Length[rowsAtD31],
                "FeasiblePoints" -> Count[Lookup[rowsAtD31, "Feasible", False], True],
                "InfeasiblePoints" -> Count[Lookup[rowsAtD31, "Feasible", False], False],
                "FeasibleIntervals" -> intervals,
                "FeasibleIntervalsText" -> ToString[InputForm[intervals]]
            |>,
            {d31, d31Values}
        ]
    ];

BuildCostMatrix[costRows_List] :=
    Module[{d21Values, d31Values, lookup, matrix},
        d21Values = Sort @ DeleteDuplicates[Lookup[costRows, "Delta21", Missing["NotAvailable"]]];
        d31Values = Sort @ DeleteDuplicates[Lookup[costRows, "Delta31", Missing["NotAvailable"]]];
        lookup = Association[
            ({Lookup[#, "Delta31"], Lookup[#, "Delta21"]} -> If[TrueQ[Lookup[#, "Feasible", False]], 1, 0]) & /@ costRows
        ];
        matrix = Join[
            {Prepend[d21Values, "Delta31\\Delta21"]},
            Table[
                Prepend[
                    (Lookup[lookup, Key[{d31, #}], ""] &) /@ d21Values,
                    d31
                ],
                {d31, d31Values}
            ]
        ];
        matrix
    ];

$ProfileConfigs = <|
    "quick" -> <|
        "Profile" -> "quick",
        "Timeout" -> 120,
        "TotalInflow" -> 258,
        "BaseI1" -> 130,
        "BaseI2" -> 128,
        "BaseU1" -> 20,
        "BaseU2" -> 100,
        "BaseU3" -> 0,
        "CostAnchorU1" -> 0,
        "FlowI1Range" -> {90, 150, 2},
        "Delta21Range" -> {40, 120, 20},
        "Delta31Range" -> {-40, 40, 20},
        "U1Range" -> {0, 140, 2},
        "U2Range" -> {-40, 140, 2},
        "U3Range" -> {-20, 80, 2}
    |>,
    "paper" -> <|
        "Profile" -> "paper",
        "Timeout" -> 120,
        "TotalInflow" -> 258,
        "BaseI1" -> 130,
        "BaseI2" -> 128,
        "BaseU1" -> 20,
        "BaseU2" -> 100,
        "BaseU3" -> 0,
        "CostAnchorU1" -> 0,
        "FlowI1Range" -> {0, 258, 1},
        "Delta21Range" -> {30, 130, 5},
        "Delta31Range" -> {-50, 50, 5},
        "U1Range" -> {0, 140, 1},
        "U2Range" -> {-40, 140, 1},
        "U3Range" -> {-20, 80, 1}
    |>,
    "dense" -> <|
        "Profile" -> "dense",
        "Timeout" -> 180,
        "TotalInflow" -> 258,
        "BaseI1" -> 130,
        "BaseI2" -> 128,
        "BaseU1" -> 20,
        "BaseU2" -> 100,
        "BaseU3" -> 0,
        "CostAnchorU1" -> 0,
        "FlowI1Range" -> {0, 258, 1},
        "Delta21Range" -> {20, 140, 2},
        "Delta31Range" -> {-60, 60, 2},
        "U1Range" -> {-20, 180, 1},
        "U2Range" -> {-60, 180, 1},
        "U3Range" -> {-40, 120, 1}
    |>
|>;

rawArgs = ParseArgs[Rest[$ScriptCommandLine]];
profile = Lookup[rawArgs, "profile", "quick"];
config = Lookup[$ProfileConfigs, profile, Lookup[$ProfileConfigs, "quick"]];
config = Join[config, <|"Profile" -> profile|>];

If[!KeyExistsQ[$ProfileConfigs, profile],
    Print["Warning: unknown profile '", profile, "'. Falling back to quick defaults with profile tag '", profile, "'."]
];

config = ApplyNumericOverride[config, rawArgs, "timeout", "Timeout"];
config = ApplyNumericOverride[config, rawArgs, "total", "TotalInflow"];
config = ApplyNumericOverride[config, rawArgs, "basei1", "BaseI1"];
config = ApplyNumericOverride[config, rawArgs, "basei2", "BaseI2"];
config = ApplyNumericOverride[config, rawArgs, "baseu1", "BaseU1"];
config = ApplyNumericOverride[config, rawArgs, "baseu2", "BaseU2"];
config = ApplyNumericOverride[config, rawArgs, "baseu3", "BaseU3"];
config = ApplyNumericOverride[config, rawArgs, "anchoru1", "CostAnchorU1"];

config = ApplyRangeOverride[config, rawArgs, "flow", "FlowI1Range"];
config = ApplyRangeOverride[config, rawArgs, "d21", "Delta21Range"];
config = ApplyRangeOverride[config, rawArgs, "d31", "Delta31Range"];
config = ApplyRangeOverride[config, rawArgs, "u1", "U1Range"];
config = ApplyRangeOverride[config, rawArgs, "u2", "U2Range"];
config = ApplyRangeOverride[config, rawArgs, "u3", "U3Range"];

Print["Configuration:"];
Print["  Profile:          ", Lookup[config, "Profile"]];
Print["  Timeout (seconds): ", Lookup[config, "Timeout"]];
Print["  Flow split:       I1 in ", Lookup[config, "FlowI1Range"], ", I2 = total - I1 with total=", Lookup[config, "TotalInflow"]];
Print["  Cost grid:        Delta21 in ", Lookup[config, "Delta21Range"], ", Delta31 in ", Lookup[config, "Delta31Range"], ", anchor U1=", Lookup[config, "CostAnchorU1"]];
Print["  U1 window:        ", Lookup[config, "U1Range"], " (I1=", Lookup[config, "BaseI1"], ", I2=", Lookup[config, "BaseI2"], ", U2=", Lookup[config, "BaseU2"], ", U3=", Lookup[config, "BaseU3"], ")"];
Print["  U2 window:        ", Lookup[config, "U2Range"], " (I1=", Lookup[config, "BaseI1"], ", I2=", Lookup[config, "BaseI2"], ", U1=", Lookup[config, "BaseU1"], ", U3=", Lookup[config, "BaseU3"], ")"];
Print["  U3 window:        ", Lookup[config, "U3Range"], " (I1=", Lookup[config, "BaseI1"], ", I2=", Lookup[config, "BaseI2"], ", U1=", Lookup[config, "BaseU1"], ", U2=", Lookup[config, "BaseU2"], ")"];
Print[];

flowI1Values = Select[
    RangeSpecToList[Lookup[config, "FlowI1Range"]],
    0 <= # <= Lookup[config, "TotalInflow"] &
];
d21Values = RangeSpecToList[Lookup[config, "Delta21Range"]];
d31Values = RangeSpecToList[Lookup[config, "Delta31Range"]];
u1Values = RangeSpecToList[Lookup[config, "U1Range"]];
u2Values = RangeSpecToList[Lookup[config, "U2Range"]];
u3Values = RangeSpecToList[Lookup[config, "U3Range"]];

flowScenarios =
    (<|
        "Slice" -> "FlowSplit",
        "SweepVariable" -> "I1",
        "SweepValue" -> #,
        "I1" -> #,
        "I2" -> Lookup[config, "TotalInflow"] - #,
        "U1" -> Lookup[config, "BaseU1"],
        "U2" -> Lookup[config, "BaseU2"],
        "U3" -> Lookup[config, "BaseU3"]
    |>&) /@ flowI1Values;

costScenarios =
    Flatten[
        Table[
            <|
                "Slice" -> "CostDifferenceGrid",
                "SweepVariable" -> "Delta21Delta31",
                "SweepValue" -> ToString[InputForm[{d21, d31}]],
                "I1" -> Lookup[config, "BaseI1"],
                "I2" -> Lookup[config, "BaseI2"],
                "U1" -> Lookup[config, "CostAnchorU1"],
                "U2" -> Lookup[config, "CostAnchorU1"] + d21,
                "U3" -> Lookup[config, "CostAnchorU1"] + d31
            |>,
            {d31, d31Values},
            {d21, d21Values}
        ],
        1
    ];

u1Scenarios =
    (<|
        "Slice" -> "U1Window",
        "SweepVariable" -> "U1",
        "SweepValue" -> #,
        "I1" -> Lookup[config, "BaseI1"],
        "I2" -> Lookup[config, "BaseI2"],
        "U1" -> #,
        "U2" -> Lookup[config, "BaseU2"],
        "U3" -> Lookup[config, "BaseU3"]
    |>&) /@ u1Values;

u2Scenarios =
    (<|
        "Slice" -> "U2Window",
        "SweepVariable" -> "U2",
        "SweepValue" -> #,
        "I1" -> Lookup[config, "BaseI1"],
        "I2" -> Lookup[config, "BaseI2"],
        "U1" -> Lookup[config, "BaseU1"],
        "U2" -> #,
        "U3" -> Lookup[config, "BaseU3"]
    |>&) /@ u2Values;

u3Scenarios =
    (<|
        "Slice" -> "U3Window",
        "SweepVariable" -> "U3",
        "SweepValue" -> #,
        "I1" -> Lookup[config, "BaseI1"],
        "I2" -> Lookup[config, "BaseI2"],
        "U1" -> Lookup[config, "BaseU1"],
        "U2" -> Lookup[config, "BaseU2"],
        "U3" -> #
    |>&) /@ u3Values;

flowRows = EvaluateSlice["FlowSplit", flowScenarios, Lookup[config, "Timeout"]];
costRows = EvaluateSlice["CostDifferenceGrid", costScenarios, Lookup[config, "Timeout"]];
u1Rows = EvaluateSlice["U1Window", u1Scenarios, Lookup[config, "Timeout"]];
u2Rows = EvaluateSlice["U2Window", u2Scenarios, Lookup[config, "Timeout"]];
u3Rows = EvaluateSlice["U3Window", u3Scenarios, Lookup[config, "Timeout"]];

allRows = Join[flowRows, costRows, u1Rows, u2Rows, u3Rows];

boundaryRows = Join[
    {
        BoundaryRow[
            flowRows,
            "FlowSplit",
            "I1",
            <|
                "FixedTotalInflow" -> Lookup[config, "TotalInflow"],
                "FixedU1" -> Lookup[config, "BaseU1"],
                "FixedU2" -> Lookup[config, "BaseU2"],
                "FixedU3" -> Lookup[config, "BaseU3"]
            |>
        ]
    },
    CostGridBoundaryRows[costRows],
    {
        BoundaryRow[
            u1Rows,
            "U1Window",
            "U1",
            <|
                "FixedI1" -> Lookup[config, "BaseI1"],
                "FixedI2" -> Lookup[config, "BaseI2"],
                "FixedU2" -> Lookup[config, "BaseU2"],
                "FixedU3" -> Lookup[config, "BaseU3"]
            |>
        ],
        BoundaryRow[
            u2Rows,
            "U2Window",
            "U2",
            <|
                "FixedI1" -> Lookup[config, "BaseI1"],
                "FixedI2" -> Lookup[config, "BaseI2"],
                "FixedU1" -> Lookup[config, "BaseU1"],
                "FixedU3" -> Lookup[config, "BaseU3"]
            |>
        ],
        BoundaryRow[
            u3Rows,
            "U3Window",
            "U3",
            <|
                "FixedI1" -> Lookup[config, "BaseI1"],
                "FixedI2" -> Lookup[config, "BaseI2"],
                "FixedU1" -> Lookup[config, "BaseU1"],
                "FixedU2" -> Lookup[config, "BaseU2"]
            |>
        ]
    }
];

costMatrix = BuildCostMatrix[costRows];

resultsDir = FileNameJoin[{Directory[], "Results"}];
If[!DirectoryQ[resultsDir], CreateDirectory[resultsDir]];

timestamp = DateString[{"Year", "Month", "Day", "-", "Hour", "Minute", "Second"}];

pointsCSV = FileNameJoin[{resultsDir, "jamarat_phase_points_" <> timestamp <> ".csv"}];
pointsJSON = FileNameJoin[{resultsDir, "jamarat_phase_points_" <> timestamp <> ".json"}];
latestPointsCSV = FileNameJoin[{resultsDir, "jamarat_phase_points_latest.csv"}];
latestPointsJSON = FileNameJoin[{resultsDir, "jamarat_phase_points_latest.json"}];

boundariesCSV = FileNameJoin[{resultsDir, "jamarat_phase_boundaries_" <> timestamp <> ".csv"}];
boundariesJSON = FileNameJoin[{resultsDir, "jamarat_phase_boundaries_" <> timestamp <> ".json"}];
latestBoundariesCSV = FileNameJoin[{resultsDir, "jamarat_phase_boundaries_latest.csv"}];
latestBoundariesJSON = FileNameJoin[{resultsDir, "jamarat_phase_boundaries_latest.json"}];

costMatrixCSV = FileNameJoin[{resultsDir, "jamarat_phase_costdiff_matrix_" <> timestamp <> ".csv"}];
latestCostMatrixCSV = FileNameJoin[{resultsDir, "jamarat_phase_costdiff_matrix_latest.csv"}];

metadataJSON = FileNameJoin[{resultsDir, "jamarat_phase_metadata_" <> timestamp <> ".json"}];
latestMetadataJSON = FileNameJoin[{resultsDir, "jamarat_phase_metadata_latest.json"}];

ExportResultsCSV[allRows, pointsCSV];
ExportResultsJSON[allRows, pointsJSON];
ExportResultsCSV[allRows, latestPointsCSV];
ExportResultsJSON[allRows, latestPointsJSON];

ExportResultsCSV[boundaryRows, boundariesCSV];
ExportResultsJSON[boundaryRows, boundariesJSON];
ExportResultsCSV[boundaryRows, latestBoundariesCSV];
ExportResultsJSON[boundaryRows, latestBoundariesJSON];

Export[costMatrixCSV, costMatrix, "CSV"];
Export[latestCostMatrixCSV, costMatrix, "CSV"];
Print["Results exported to ", costMatrixCSV];
Print["Results exported to ", latestCostMatrixCSV];

metadata = <|
    "GeneratedAt" -> DateString[],
    "Profile" -> Lookup[config, "Profile"],
    "Configuration" -> config,
    "Counts" -> <|
        "TotalPoints" -> Length[allRows],
        "FeasiblePoints" -> Count[Lookup[allRows, "Feasible", False], True],
        "InfeasiblePoints" -> Count[Lookup[allRows, "Feasible", False], False],
        "FailedOrTimedOut" -> Count[Lookup[allRows, "ExecutionStatus", "UNKNOWN"], "FAILED" | "TIMEOUT"]
    |>,
    "SliceCounts" -> Association[
        KeyValueMap[
            #1 -> <|
                "Points" -> Length[#2],
                "Feasible" -> Count[Lookup[#2, "Feasible", False], True],
                "Infeasible" -> Count[Lookup[#2, "Feasible", False], False]
            |>&,
            GroupBy[allRows, Lookup[#, "Slice", "Unknown"] &]
        ]
    ]
|>;

Export[metadataJSON, JSONSafeValue[metadata], "JSON"];
Export[latestMetadataJSON, JSONSafeValue[metadata], "JSON"];
Print["Results exported to ", metadataJSON];
Print["Results exported to ", latestMetadataJSON];

Print[];
Print["Summary:"];
Print["  Total points:      ", metadata["Counts", "TotalPoints"]];
Print["  Feasible points:   ", metadata["Counts", "FeasiblePoints"]];
Print["  Infeasible points: ", metadata["Counts", "InfeasiblePoints"]];
Print["  Failed/timeout:    ", metadata["Counts", "FailedOrTimedOut"]];
Print[];
Print["Done."];
