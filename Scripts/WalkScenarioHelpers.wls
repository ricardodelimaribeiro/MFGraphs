(* WalkScenarioHelpers.wls -- the gate-by-gate scenario walk shared by
   Scripts/WalkScenario.wls (the CLI behind the /scenario skill) and
   MFGraphs/Tests/walk-scenario.mt (fast-suite enforcement).

   Sourced with Get[] after MFGraphs` is loaded. Defines functions only:
   no Exit, no Print, no file writes. Uses public MFGraphs symbols only, so
   an API change breaks the test instead of silently outdating the skill.

   Defines:
     walkScenario[input, opts]  - run the six gates and return an Association
     formatWalk[walk]           - render the fixed text report

   The gates, in order; the walk stops at the first one that does not pass:
     scenario  scenarioQ and validateScenario
     unknowns  symbolicUnknownsQ
     system    mfgSystemQ and criticalCongestionSystemQ
     solve     dnfReduceSystem returns a solution within the timeout
     valid     isValidSystemSolution is True
     exact     exactSolutionReport proves soundness

   Gate statuses: PASS, FAIL, NOSOLUTION, TIMEOUT, SKIP (not reached).
   Results: "OK", "GateFailed", "NoSolution", "Timeout". *)

$walkGateNames = {"scenario", "unknowns", "system", "solve", "valid", "exact"};

(* The reference solver recurses deeply on large grids. At the default limit of
   1024 the kernel answers TerminatedEvaluation[RecursionLimit], which no
   Catch, CheckAbort or TimeConstrained can intercept: a script simply stops
   and wolframscript still exits 0. The walk therefore runs with the same
   raised limits as MFGraphs/Tests/example-coverage.mt. *)
$walkRecursionLimit = 16384;

Options[walkScenario] = {
    "Timeout"      -> 60,     (* seconds allowed for the solve gate *)
    "Completeness" -> None    (* seconds for the completeness search; None skips it *)
};

walkClip[text_String] :=
    If[StringLength[text] > 400, StringTake[text, 400] <> " ...", text];

walkFailureDetail[failure_Failure] :=
    Module[{tag, data, message, rest},
        tag = First[failure];
        data = If[Length[failure] >= 2 && AssociationQ[failure[[2]]], failure[[2]], <||>];
        message = Lookup[data, "Message", failure["Message"]];
        rest = Select[
            KeyDrop[data, {"Message", "MessageTemplate", "MessageParameters"}],
            # =!= {} && # =!= <||> &
        ];
        walkClip[
            ToString[tag] <> ": " <>
            If[StringQ[message], message, ToString[message, InputForm]] <>
            If[Length[rest] > 0, " " <> ToString[Normal[rest], InputForm], ""]
        ]
    ];
walkFailureDetail[other_] :=
    walkClip["unexpected value " <> ToString[Short[other, 2], InputForm]];

walkSeconds[seconds_] :=
    ToString[NumberForm[N[seconds], {Infinity, 2}, NumberPadding -> {"", "0"}]] <> " s";

walkNoSolutionQ[solution_] :=
    solution === {} ||
    (AssociationQ[solution] && Lookup[solution, "Residual", True] === False);

walkSolutionShape[solution_] :=
    Which[
        ListQ[solution],
            ToString[Length[solution]] <> " rules, no residual",
        Lookup[solution, "Residual", True] === True,
            ToString[Length[solution["Rules"]]] <> " rules, no residual",
        True,
            ToString[Length[solution["Rules"]]] <> " rules + residual"
    ];

walkScenario[input_, OptionsPattern[]] :=
    Module[{gates = {}, addGate, finishWalk, timeout, completeness,
            scen = Missing["NotReached"], unknownsBundle = Missing["NotReached"],
            sys = Missing["NotReached"], solution = Missing["NotReached"],
            exact = Missing["NotReached"], model, seconds, report, valid, exactDetail},
        timeout = OptionValue["Timeout"];
        completeness = OptionValue["Completeness"];

        addGate[name_, status_, detail_] :=
            AppendTo[gates, <|"Gate" -> name, "Status" -> status, "Detail" -> detail|>];
        finishWalk[result_] := (
            Scan[
                If[!MemberQ[#["Gate"] & /@ gates, #], addGate[#, "SKIP", ""]] &,
                $walkGateNames
            ];
            Throw[
                <|"Gates" -> gates, "Result" -> result, "Scenario" -> scen,
                  "System" -> sys, "Solution" -> solution, "Exact" -> exact|>,
                $walkFinished
            ]
        );

        Catch[
          Block[{$RecursionLimit = $walkRecursionLimit, $IterationLimit = $walkRecursionLimit},
            (* Gate 1: scenario *)
            scen = Which[
                scenarioQ[input], input,
                AssociationQ[input], makeScenario[input],
                True, input
            ];
            If[!scenarioQ[scen],
                addGate["scenario", "FAIL", walkFailureDetail[scen]];
                finishWalk["GateFailed"]
            ];
            If[!TrueQ[validateScenario[scen]],
                addGate["scenario", "FAIL",
                    "validateScenario returned False; its message names the first problem"];
                finishWalk["GateFailed"]
            ];
            model = scenarioData[scen, "Model"];
            addGate["scenario", "PASS",
                ToString[Length[model["Vertices"]]] <> " vertices, entries " <>
                ToString[model["Entries"], InputForm] <> ", exits " <>
                ToString[model["Exits"], InputForm]];

            (* Gate 2: unknowns *)
            unknownsBundle = makeSymbolicUnknowns[scen];
            If[!symbolicUnknownsQ[unknownsBundle],
                addGate["unknowns", "FAIL", walkFailureDetail[unknownsBundle]];
                finishWalk["GateFailed"]
            ];
            addGate["unknowns", "PASS",
                StringRiffle[
                    # <> " " <> ToString[Length[symbolicUnknownsData[unknownsBundle, #]]] & /@
                        {"Js", "Jts", "Us"},
                    ", "
                ]];

            (* Gate 3: system *)
            sys = makeSystem[scen, unknownsBundle];
            If[!mfgSystemQ[sys],
                addGate["system", "FAIL", walkFailureDetail[sys]];
                finishWalk["GateFailed"]
            ];
            If[!TrueQ[criticalCongestionSystemQ[sys]],
                addGate["system", "FAIL",
                    "not critical congestion: the solvers accept Alpha = 1 on every edge only"];
                finishWalk["GateFailed"]
            ];
            addGate["system", "PASS", "critical congestion (Alpha = 1 on every edge)"];

            (* Gate 4: solve. The timeout marker is tested here, before any
               report function sees the result: the package does not recognise
               a caller's own marker symbol. *)
            {seconds, solution} = AbsoluteTiming[
                TimeConstrained[dnfReduceSystem[sys], timeout, $walkTimedOut]
            ];
            Which[
                solution === $walkTimedOut,
                    addGate["solve", "TIMEOUT",
                        "dnfReduceSystem exceeded " <> ToString[timeout] <>
                        " s; this says nothing about feasibility"];
                    finishWalk["Timeout"],
                FailureQ[solution],
                    addGate["solve", "FAIL", walkFailureDetail[solution]];
                    finishWalk["GateFailed"],
                !(ListQ[solution] || (AssociationQ[solution] && KeyExistsQ[solution, "Rules"])),
                    addGate["solve", "FAIL", walkFailureDetail[solution]];
                    finishWalk["GateFailed"],
                walkNoSolutionQ[solution],
                    addGate["solve", "NOSOLUTION",
                        "dnfReduceSystem " <> walkSeconds[seconds] <>
                        ": the residual reduced to False"];
                    finishWalk["NoSolution"]
            ];
            report = solutionReport[sys, solution];
            addGate["solve", "PASS",
                "dnfReduceSystem " <> walkSeconds[seconds] <> ", " <>
                walkSolutionShape[solution] <> ", ResultKind " <>
                ToString[Lookup[report, "ResultKind", "Unknown"]]];

            (* Gate 5: valid *)
            valid = isValidSystemSolution[sys, solution];
            If[valid =!= True,
                addGate["valid", "FAIL",
                    "isValidSystemSolution returned " <> ToString[valid, InputForm]];
                finishWalk["GateFailed"]
            ];
            addGate["valid", "PASS",
                "isValidSystemSolution True, Kirchhoff residual " <>
                ToString[Lookup[report, "KirchhoffResidual", Missing["NotAvailable"]], InputForm]];

            (* Gate 6: exact soundness; completeness only on request *)
            exact = If[NumericQ[completeness] && completeness > 0,
                exactSolutionReport[sys, solution, "TimeLimit" -> completeness],
                exactSolutionReport[sys, solution, "CheckCompleteness" -> False]
            ];
            If[!AssociationQ[exact],
                addGate["exact", "FAIL", walkFailureDetail[exact]];
                finishWalk["GateFailed"]
            ];
            exactDetail = StringRiffle[
                # <> " " <> ToString[Lookup[exact, #, Missing["KeyAbsent", #]]] & /@
                    {"ResultType", "Soundness", "Completeness"},
                ", "
            ];
            If[exact["Soundness"] =!= "Proved",
                addGate["exact", "FAIL", exactDetail];
                finishWalk["GateFailed"]
            ];
            addGate["exact", "PASS", exactDetail];
            finishWalk["OK"]
          ]
          ,
          $walkFinished
        ]
    ];

formatWalk[walk_Association] :=
    StringRiffle[
        Join[
            MapIndexed[
                Function[{gate, position},
                    StringTrim[
                        "[" <> ToString[First[position]] <> "] " <>
                        StringPadRight[gate["Gate"], 9] <> " " <>
                        StringPadRight[gate["Status"], 11] <> gate["Detail"]
                    ]
                ],
                walk["Gates"]
            ],
            {"RESULT        " <> walk["Result"]}
        ],
        "\n"
    ];
