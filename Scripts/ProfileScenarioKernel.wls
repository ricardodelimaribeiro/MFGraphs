#!/usr/bin/env wolframscript
(* Usage: wolframscript -file Scripts/ProfileScenarioKernel.wls [--solver dnf|optimizeddnf|activeset|reduce|boolean|findinstance] [--timeout s] [--case name] *)
(* Non-mutating staged profiler for scenario -> symbolic unknowns -> system -> solver. *)

(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

$RepoRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
PrependTo[$Path, $RepoRoot];
Needs["MFGraphs`"];
$MFGraphsVerbose = False;
$HistoryLength = 0;

$TimeoutSecs = 60;
$SelectedCase = All;
$SolverName = "dnf";

Module[{args = Rest[$ScriptCommandLine], i = 1, arg},
  While[i <= Length[args],
    arg = args[[i]];
    Switch[arg,
      "--timeout",
        i++; If[i <= Length[args],
          Module[{v = Quiet @ Check[ToExpression[args[[i]]], $Failed]},
            If[NumericQ[v] && v > 0, $TimeoutSecs = v]]],
      "--case",
        i++; If[i <= Length[args], $SelectedCase = args[[i]]],
      "--solver",
        i++; If[i <= Length[args], $SolverName = ToLowerCase[args[[i]]]],
      _, Null
    ];
    i++;
  ]
];

$SolverNames = {"dnf", "optimizeddnf", "activeset", "reduce", "boolean", "findinstance"};

If[!MemberQ[$SolverNames, $SolverName],
  Print["Unknown solver: ", $SolverName];
  Print["Available solvers: ", StringRiffle[$SolverNames, ", "]];
  Exit[2]
];

RunSelectedSolver[sys_, constraints_, allVars_, rulesAcc_] :=
  Switch[$SolverName,
    "dnf", solversTools`dnfReduceSystem[sys],
    "optimizeddnf", solversTools`optimizedDNFReduceSystem[sys],
    "activeset", solversTools`activeSetReduceSystem[sys],
    "reduce",
      With[{reduced = Reduce[constraints, allVars, Reals]},
        solversTools`Private`attachAccumulatedRules[
          solversTools`Private`parseReduceResult[reduced, allVars],
          rulesAcc
        ]
      ],
    "boolean", solversTools`booleanReduceSystem[sys],
    "findinstance", solversTools`findInstanceSystem[sys]
  ];

Fmt[x_?NumericQ] := ToString[NumberForm[N[x], {Infinity, 3}, NumberPadding -> {"", "0"}]];
Fmt[x_] := ToString[x, InputForm];
Ms[t_?NumericQ] := N[1000 t];
AndLen[expr_] := If[Head[expr] === And, Length[List @@ expr], If[TrueQ[expr], 0, 1]];

ResultKind[sol_] :=
  With[{kind = solversTools`Private`solutionResultKind[sol]},
    If[kind =!= "Unknown", kind, ToString[Head[sol]]]
  ];

ModelEdgeCount[model_Association] :=
  Module[{graph = Lookup[model, "Graph", Missing["KeyAbsent", "Graph"]],
          adjacency = Lookup[model, "Adjacency", Missing["KeyAbsent", "Adjacency"]]},
    Which[
      MatchQ[graph, _Graph], EdgeCount[graph],
      MatrixQ[adjacency], Total[Flatten[Unitize[adjacency]]],
      True, Missing["Unavailable"]
    ]
  ];

$Cases = {
  {"chain-2v", Function[gridScenario[{2}, {{1, 10}}, {{2, 0}}]]},
  {"chain-3v-1exit", Function[gridScenario[{3}, {{1, 10}}, {{3, 0}}]]},
  {"chain-3v-2exit", Function[gridScenario[{3}, {{1, 120.0}}, {{2, 0.0}, {3, 10.0}}]]},
  {"chain-3-midentry", Function[gridScenario[{3}, {{2, 100.0}}, {{1, 0.0}, {3, 10.0}}]]},
  {"chain-5v-1exit", Function[gridScenario[{5}, {{1, 10}}, {{5, 0}}]]},
  {"grid-2x3", Function[gridScenario[{2, 3}, {{1, 100}}, {{6, 0}}]]},
  {"grid-3x2", Function[gridScenario[{3, 2}, {{1, 100}}, {{6, 0}}]]},
  {"example-12", Function[getExampleScenario[12, {{1, 100.0}}, {{4, 0.0}}]]}
};
$CaseNames = $Cases[[All, 1]];

If[$SelectedCase =!= All,
  $Cases = Select[$Cases, #[[1]] === $SelectedCase &];
  If[$Cases === {},
    Print["No case named '", $SelectedCase, "'. Available: ",
      StringRiffle[$CaseNames, ", "]];
    Exit[2]
  ]
];

ProfileCase[{name_String, buildFn_}] :=
  Module[{mem0, mem1, mem2, mem3, mem4, tScenario, s, tUnknowns, unk,
          tSystem, sys, topology, tInputs, inputs, constraints, allVars,
          rulesAcc, tSolve, sol, kind, valid, status, model, vertices,
          edges, auxPairs, auxTriples, halfPairs, js, jts, us, eqCount,
          switchIneqCount, altFlowCount, altTransitionCount, altOptCount,
          scenarioBytes, unknownBytes, systemBytes, varDiag, transitionStatus,
          transitionRules, transitionResiduals},
    Print["CASE ", name];

    mem0 = MemoryInUse[];
    {tScenario, s} = AbsoluteTiming[buildFn[]];
    mem1 = MemoryInUse[];

    {tUnknowns, unk} = AbsoluteTiming[makeSymbolicUnknowns[s]];
    mem2 = MemoryInUse[];

    {tSystem, sys} = AbsoluteTiming[makeSystem[s, unk]];
    mem3 = MemoryInUse[];

    model = scenarioData[s, "Model"];
    topology = scenarioData[s, "Topology"];
    vertices = Length[model["Vertices"]];
    edges = ModelEdgeCount[model];
    auxPairs = Length[Lookup[topology, "AuxPairs", {}]];
    auxTriples = Length[Lookup[topology, "AuxTriples", {}]];
    halfPairs = Length[Lookup[topology, "HalfPairs", {}]];
    js = Length[symbolicUnknownsData[unk, "Js"]];
    jts = Length[symbolicUnknownsData[unk, "Jts"]];
    us = Length[symbolicUnknownsData[unk, "Us"]];
    eqCount = Length[systemData[sys, "EqEntryIn"]] +
      AndLen[systemData[sys, "EqBalanceSplittingFlows"]] +
      AndLen[systemData[sys, "EqBalanceGatheringFlows"]] +
      AndLen[systemData[sys, "EqGeneral"]];
    switchIneqCount = AndLen[systemData[sys, "IneqSwitchingByVertex"]];
    altFlowCount = AndLen[systemData[sys, "AltFlows"]];
    altTransitionCount = AndLen[systemData[sys, "AltTransitionFlows"]];
    altOptCount = AndLen[systemData[sys, "AltOptCond"]];
    scenarioBytes = ByteCount[s];
    unknownBytes = ByteCount[unk];
    systemBytes = ByteCount[sys];

    {tInputs, inputs} = AbsoluteTiming[solversTools`Private`buildSolverInputs[sys]];
    {constraints, allVars, rulesAcc} = inputs;

    {tSolve, sol} = AbsoluteTiming[TimeConstrained[RunSelectedSolver[sys, constraints, allVars, rulesAcc], $TimeoutSecs, $TimedOut]];
    mem4 = MemoryInUse[];
    kind = ResultKind[sol];
    status = If[sol === $TimedOut, "TIMEOUT", "OK"];
    varDiag = If[sol === $TimedOut, <||>, solversTools`Private`solutionVariableDiagnostics[sys, sol]];
    transitionStatus = Lookup[varDiag, "TransitionFlowStatus", "TimedOut"];
    transitionRules = Lookup[varDiag, "TransitionFlowRuleCount", "TimedOut"];
    transitionResiduals = Lookup[varDiag, "TransitionFlowResidualCount", "TimedOut"];
    valid = If[sol === $TimedOut, "TimedOut", ToString[isValidSystemSolution[sys, sol]]];

    Print[StringRiffle[
      {
        name,
        Fmt[Ms[tScenario]], Fmt[Ms[tUnknowns]], Fmt[Ms[tSystem]], Fmt[Ms[tInputs]], Fmt[Ms[tSolve]],
        status, kind, transitionStatus, ToString[transitionRules], ToString[transitionResiduals], valid,
        ToString[vertices], ToString[edges], ToString[halfPairs], ToString[auxPairs], ToString[auxTriples],
        ToString[js], ToString[jts], ToString[us], ToString[Length[allVars]], ToString[Length[rulesAcc]],
        ToString[eqCount], ToString[switchIneqCount], ToString[altFlowCount], ToString[altTransitionCount], ToString[altOptCount],
        ToString[scenarioBytes], ToString[unknownBytes], ToString[systemBytes],
        ToString[mem1 - mem0], ToString[mem2 - mem1], ToString[mem3 - mem2], ToString[mem4 - mem3]
      },
      "\t"
    ]];
  ];

Print["MFGraphs staged scenario-kernel profile"];
Print["Date\t", DateString[]];
Print["Solver\t", $SolverName];
Print["TimeoutSeconds\t", $TimeoutSecs];
If[$SelectedCase =!= All, Print["CaseFilter\t", $SelectedCase]];
Print["Header\tCase\tScenarioMs\tUnknownsMs\tSystemMs\tSolverInputMs\tSolveMs\tStatus\tKind\tTransitionFlowStatus\tTransitionFlowRules\tTransitionFlowResiduals\tValid\tVertices\tEdges\tHalfPairs\tAuxPairs\tAuxTriples\tJs\tJts\tUs\tSolverVars\tPreRules\tEqCount\tSwitchIneqs\tAltFlows\tAltTransitionFlows\tAltOptCond\tScenarioBytes\tUnknownBytes\tSystemBytes\tScenarioMemDelta\tUnknownsMemDelta\tSystemMemDelta\tSolveMemDelta"];
Scan[ProfileCase, $Cases];

(* :!CodeAnalysis::EndBlock:: *)
